#!/usr/bin/env python3
"""Resumable nested-prefix benchmark; source datasets are only read."""
import argparse
import concurrent.futures
import contextlib
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import signal
import statistics
import subprocess
import time


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_bytes(path):
    total = 0
    for directory, _, files in os.walk(path):
        for name in files:
            try:
                total += (Path(directory) / name).stat().st_size
            except FileNotFoundError:
                pass  # Construction deletes its own temporary buckets.
    return total


def available_memory():
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) * 1024
    raise RuntimeError("Cannot read available memory")


def memory_headroom(available=None, arcstats=Path("/proc/spl/kstat/zfs/arcstats")):
    """Discount clean ZFS data pages; ARC is outside Linux's page cache.

    Never credit ghost lists, dirty/referenced data, metadata, or the ARC
    minimum. Credit only half the remaining evictable data as a margin for
    reclaim latency. No system/cache configuration is changed.
    """
    available = available_memory() if available is None else available
    result = {"mem_available_bytes": available, "zfs_arc_bytes": 0,
              "zfs_evictable_data_bytes": 0, "zfs_arc_min_bytes": 0,
              "zfs_credit_bytes": 0, "effective_available_bytes": available}
    try:
        stats = {}
        for line in arcstats.read_text().splitlines():
            fields = line.split()
            if len(fields) == 3:
                try:
                    stats[fields[0]] = int(fields[2])
                except ValueError:
                    pass
        size, minimum = stats["size"], stats["c_min"]
        clean = stats["mru_evictable_data"] + stats["mfu_evictable_data"]
        if min(size, minimum, clean) < 0:
            return result
        credit = min(clean, max(0, size - minimum)) // 2
        result.update(zfs_arc_bytes=size, zfs_evictable_data_bytes=clean,
                      zfs_arc_min_bytes=minimum, zfs_credit_bytes=credit,
                      effective_available_bytes=available + credit)
    except (OSError, KeyError):
        pass  # Preserve the Linux-only guard if ARC statistics are unavailable.
    return result


def memory_emergency(headroom):
    return headroom["mem_available_bytes"] < 8 * 1024**3


def process_resources(pid):
    """Own process tree only; Linux counters are cumulative since process start."""
    pending = [pid]
    rss = read_bytes = write_bytes = 0
    while pending:
        current = pending.pop()
        directory = Path(f"/proc/{current}")
        try:
            children = directory / f"task/{current}/children"
            pending.extend(map(int, children.read_text().split()))
            for line in (directory / "status").read_text().splitlines():
                if line.startswith("VmRSS:"): rss += int(line.split()[1]) * 1024
            counters = dict(line.split(": ") for line in (directory / "io").read_text().splitlines())
            read_bytes += int(counters["read_bytes"])
            write_bytes += int(counters["write_bytes"])
        except (OSError, ValueError, KeyError):
            pass  # A measured child may finish between counter reads.
    return {"rss_bytes": rss, "read_bytes": read_bytes, "write_bytes": write_bytes}


def measure(command, log, index=None):
    timing = log.with_suffix(".time")
    start = time.monotonic()
    peak_disk = 0
    samples = []
    resource_limit = None
    working = index if index else log.parent.parent
    working.mkdir(parents=True, exist_ok=True)
    with log.open("w") as output:
        process = subprocess.Popen(["/usr/bin/time", "-f", "%e %U %S %M", "-o",
                                    str(timing), *map(str, command)],
                                   stdout=output, stderr=subprocess.STDOUT,
                                   start_new_session=True, cwd=working)
        last_disk = 0
        low_memory = 0
        while process.poll() is None:
            now = time.monotonic()
            headroom = memory_headroom()
            if index and now - last_disk >= 20:
                disk = tree_bytes(index)
                peak_disk = max(peak_disk, disk)
                samples.append({"elapsed_seconds": now - start, "disk_bytes": disk,
                                **process_resources(process.pid), **headroom})
                last_disk = now
            low_memory = low_memory + 1 if headroom["effective_available_bytes"] < 24 * 1024**3 else 0
            filesystem = os.statvfs(log.parent)
            if low_memory >= 3:
                resource_limit = "Server RAM including discounted clean ZFS ARC below 24 GiB for 15 seconds"
            if memory_emergency(headroom):
                resource_limit = "Linux MemAvailable below the 8 GiB emergency floor"
            if filesystem.f_bavail * filesystem.f_frsize < 100 * 1024**3:
                resource_limit = "Dedicated filesystem has less than 100 GiB free"
            if resource_limit:
                # Let GNU time reap its child and retain resource measurements.
                children = Path(f"/proc/{process.pid}/task/{process.pid}/children")
                for pid in children.read_text().split() if children.exists() else []:
                    try:
                        os.kill(int(pid), signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                try:
                    process.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                break
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
        code = process.wait()
    result = {"command": list(map(str, command)), "exit_code": code,
              "wall_seconds": time.monotonic() - start,
              "peak_construction_disk_bytes_sampled": peak_disk,
              "disk_sampling_interval_seconds": 20, "resource_limit": resource_limit,
              "memory_guard_policy": "24 GiB with half clean ARC credit, capped above c_min; 8 GiB Linux emergency floor"}
    if index:
        result["temporary_roots_measured"] = [str(index)]
        result["construction_resource_samples"] = samples
    if timing.exists():
        lines = timing.read_text().splitlines()
        fields = lines[-1].split() if lines else []
        if len(fields) == 4:
            result.update(elapsed_seconds=float(fields[0]), user_seconds=float(fields[1]),
                          system_seconds=float(fields[2]), peak_rss_kib=int(fields[3]))
    return result


def manifest(base, root, maximum):
    path = base / "reports/datasets.json"
    if path.exists():
        selected = json.loads(path.read_text())
        assert len(selected) >= maximum
        return selected
    rng = random.Random(381592)
    folders = sorted(p for p in root.iterdir() if p.is_dir() and p.name.isdigit())
    rng.shuffle(folders)
    choices = []
    for folder in folders:
        files = sorted(folder.glob("*.ru.fa.zstd"))
        rng.shuffle(files)
        choices.append(files)
    selected = []
    round_number = 0
    names = set()
    while len(selected) < maximum:
        progress = False
        for files in choices:
            if round_number < len(files):
                original = files[round_number]
                if original.name in names:
                    continue
                names.add(original.name)
                destination = base / "inputs" / original.parent.name / original.name.removesuffix(".zstd")
                selected.append({"source": str(original), "input": str(destination),
                                 "compressed_bytes": original.stat().st_size})
                progress = True
                if len(selected) == maximum:
                    break
        assert progress, "Not enough distinct datasets"
        round_number += 1
    save(path, selected)
    return selected


def stage(entry):
    destination = Path(entry["input"])
    if destination.exists():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".partial")
    start = time.monotonic()
    with temporary.open("wb") as output:
        subprocess.run(["zstd", "-q", "-dc", entry["source"]], stdout=output, check=True)
    temporary.replace(destination)
    save(destination.with_suffix(".staging.json"),
         {**entry, "plain_bytes": destination.stat().st_size,
          "decompression_seconds": time.monotonic() - start})


def records(path):
    with path.open("rb") as stream:
        sequence = []
        for line in stream:
            if line.startswith(b">"):
                if sequence:
                    yield b"".join(sequence)
                sequence = []
            else:
                sequence.append(line.strip())
        if sequence:
            yield b"".join(sequence)


def queries(base, entries):
    path = base / "inputs/queries.fa"
    if path.exists():
        return path
    rng = random.Random(381593)
    pool = []
    totals = []
    for dataset, entry in enumerate(entries[:16]):
        reservoir = []
        seen = 0
        for sequence in records(Path(entry["input"])):
            sequence = sequence.upper()
            if len(sequence) < 31 or sequence.strip(b"ACGT"):
                continue
            seen += 1
            slot = rng.randrange(seen)
            if len(reservoir) < 64:
                reservoir.append(sequence[:2000])
            elif slot < 64:
                reservoir[slot] = sequence[:2000]
        totals.append(seen)
        pool.extend((f"dataset{dataset}_{i}", s) for i, s in enumerate(reservoir))
    original = list(pool)
    table = bytes.maketrans(b"ACGT", b"TGCA")
    pool.extend((name + "_rc", sequence.translate(table)[::-1]) for name, sequence in original)
    for i in range(256):
        pool.append((f"random{i}", bytes(rng.choices(b"ACGT", k=150))))
    with path.open("wb") as output:
        for name, sequence in pool:
            output.write(b">" + name.encode() + b"\n" + sequence + b"\n")
    save(base / "reports/queries.json", {"seed": 381593, "records": len(pool),
         "kmers": sum(len(s) - 30 for _, s in pool), "sha256": sha(path),
         "valid_unitig_records_in_first_16_datasets": totals})
    return path


def bulk_queries(base, query):
    path = base / "inputs/queries-throughput.fa"
    if not path.exists():
        blocks = query.read_bytes().split(b">")[1:]
        with path.open("wb") as output:
            for repeat in range(64):
                for block in blocks:
                    header, sequence = block.split(b"\n", 1)
                    output.write(b">r" + str(repeat).encode() + b"_" + header + b"\n" + sequence)
        original = json.loads((base / "reports/queries.json").read_text())
        save(base / "reports/queries-throughput.json", {"records": original["records"] * 64,
             "kmers": original["kmers"] * 64, "sha256": sha(path), "repetitions": 64})
    return path


def measure_bulk(base, kind, count, threads, query, result):
    output = base / f"reports/{kind}-{count}-throughput.tsv"
    measurements = []
    for run in range(3):
        measurement = measure([base / f"bin/reindeer-{kind}", "--query", "-l", base / f"indexes/{kind}-{count}",
                               "-q", query, "-o", output, "-t", threads, "-P", 0],
                              base / f"logs/{kind}-{count}-throughput-{run}.log")
        if measurement["exit_code"] == 0:
            measurement["output_sha256"] = sha(output)
        measurements.append(measurement)
    result["bulk_queries"] = measurements
    if all(m["exit_code"] == 0 for m in measurements):
        result["bulk_query_seconds_median"] = statistics.median(m["elapsed_seconds"] for m in measurements)
        result["bulk_outputs_repeatable"] = len({m["output_sha256"] for m in measurements}) == 1
        reference = base / f"reports/cpp-{count}-throughput.tsv"
        result["bulk_outputs_match_cpp"] = sha(output) == sha(reference) if reference.exists() else None


def validate(base, count, kinds):
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, max(soft, count + 64)), hard))
    report = {"datasets": count, "oracle": "Independent canonical k=31 FASTA scan; maximum input abundance", "implementations": {}}
    for kind in kinds:
        output = base / f"reports/{kind}-{count}-query.tsv"
        if not output.exists():
            continue
        checked = 0
        mismatches = 0
        examples = []
        with contextlib.ExitStack() as stack:
            references = [iter(csv.reader(stack.enter_context((base / f"reports/oracle/{color}.tsv").open()), delimiter="\t"))
                          for color in range(count)]
            rows = csv.reader(stack.enter_context(output.open()), delimiter="\t")
            header = next(rows)
            assert len(header) == count + 1, (kind, "Incorrect dataset columns")
            for row in rows:
                assert len(row) == count + 1
                for color, reference in enumerate(references):
                    expected = next(reference)
                    checked += 1
                    if expected != [row[0], row[color + 1]]:
                        mismatches += 1
                        if len(examples) < 20:
                            examples.append({"query": row[0], "dataset": header[color + 1],
                                             "expected": expected[1], "actual": row[color + 1]})
            assert all(next(reference, None) is None for reference in references)
        report["implementations"][kind] = {"query_dataset_pairs_checked": checked,
            "mismatches": mismatches, "examples": examples, "output_sha256": sha(output)}
    save(base / f"reports/correctness-{count}.json", report)
    for kind in report["implementations"]:
        save(base / f"reports/correctness-{count}-{kind}.json", report)
    print(json.dumps(report), flush=True)
    if any(result["mismatches"] for result in report["implementations"].values()):
        raise RuntimeError(f"Independent abundance oracle mismatch at {count} datasets; see correctness report")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=Path("/DATA/alimasse"))
    parser.add_argument("--maximum", type=int, default=1024)
    parser.add_argument("--minimum", type=int, default=2)
    parser.add_argument("--threads", type=int, default=64)
    parser.add_argument("--kinds", nargs="+", default=["cpp", "rust-baseline"])
    parser.add_argument("--index-flags", nargs=argparse.REMAINDER, default=[])
    args = parser.parse_args()
    if args.maximum < 2 or args.maximum & (args.maximum - 1):
        parser.error("--maximum must be a power of two, at least 2")
    if args.minimum < 2 or args.minimum & (args.minimum - 1) or args.minimum > args.maximum:
        parser.error("--minimum must be a power of two between 2 and --maximum")
    if args.threads < 1:
        parser.error("--threads must be positive")
    base = args.base.resolve()
    for directory in ["bin", "inputs", "indexes", "logs", "reports"]:
        (base / directory).mkdir(parents=True, exist_ok=True)
    for kind in args.kinds:
        if not (base / f"bin/reindeer-{kind}").is_file():
            parser.error(f"Missing benchmark binary: {base / f'bin/reindeer-{kind}'}")
    run_config = {"minimum": args.minimum, "maximum": args.maximum,
                  "threads": args.threads, "kinds": args.kinds}
    config_path = base / ("reports/run-config-" + "_".join(args.kinds) + ".json")
    if config_path.exists() and json.loads(config_path.read_text()) != run_config:
        raise RuntimeError("Benchmark range or configuration changed; use a new variant name")
    save(config_path, run_config)
    entries = manifest(base, args.root, args.maximum)
    variants_path = base / "reports/variant-flags.json"
    variants = json.loads(variants_path.read_text()) if variants_path.exists() else {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as workers:
        list(workers.map(stage, entries[:16]))
    query = queries(base, entries)
    bulk = bulk_queries(base, query)
    results_path = base / ("reports/results-" + "_".join(args.kinds) + ".json")
    results = json.loads(results_path.read_text()) if results_path.exists() else []
    count = args.minimum
    while count <= args.maximum:
        while (base / "reports/benchmark-hold").exists():
            time.sleep(5)
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as workers:
            list(workers.map(stage, entries[:count]))
        fof = base / f"inputs/datasets-{count}.fof"
        fof.write_text("".join(entry["input"] + "\n" for entry in entries[:count]))
        for kind in args.kinds:
            previous = next((r for r in results if r["kind"] == kind and r["datasets"] == count), None)
            if previous is not None:
                if previous["binary_sha256"] != sha(base / f"bin/reindeer-{kind}") or previous["threads"] != args.threads:
                    raise RuntimeError(f"Cannot resume {kind}-{count}: binary or thread count changed; use a new variant name")
                if previous["build"]["exit_code"] != 0:
                    raise RuntimeError(f"Cannot resume failed {kind}-{count}; preserve the attempt and use a repaired variant")
                previous_flags = previous["build"]["command"][previous["build"]["command"].index("-t") + 2:]
                if previous_flags != args.index_flags + variants.get(kind, []):
                    raise RuntimeError(f"Cannot resume {kind}-{count}: construction flags changed; use a new variant name")
                if previous["build"]["exit_code"] == 0 and "bulk_queries" not in previous:
                    measure_bulk(base, kind, count, args.threads, bulk, previous)
                    save(results_path, results)
                    print(json.dumps(previous), flush=True)
                continue
            binary = base / f"bin/reindeer-{kind}"
            index = base / f"indexes/{kind}-{count}"
            result = {"kind": kind, "datasets": count, "threads": args.threads,
                      "binary_sha256": sha(binary),
                      "input_plain_bytes": sum(Path(e["input"]).stat().st_size for e in entries[:count]),
                      "input_compressed_bytes": sum(e["compressed_bytes"] for e in entries[:count])}
            result["build"] = measure([binary, "--index", "-f", fof, "-o", index,
                                        "-t", args.threads, *args.index_flags, *variants.get(kind, [])], base / f"logs/{kind}-{count}-build.log", index)
            if result["build"]["exit_code"] == 0:
                components = {str(p.relative_to(index)): p.stat().st_size
                              for p in index.rglob("*") if p.is_file()}
                result["index_components_bytes"] = components
                result["index_bytes"] = sum(components.values())
                result["queries"] = []
                for run in range(3):
                    output = base / f"reports/{kind}-{count}-query.tsv"
                    measurement = measure([binary, "--query", "-l", index, "-q", query,
                                           "-o", output, "-t", args.threads, "-P", 0],
                                          base / f"logs/{kind}-{count}-query-{run}.log")
                    if measurement["exit_code"] == 0:
                        measurement["output_sha256"] = sha(output)
                    result["queries"].append(measurement)
                if all(q["exit_code"] == 0 for q in result["queries"]):
                    result["outputs_repeatable"] = len({q["output_sha256"] for q in result["queries"]}) == 1
                    result["query_seconds_median"] = statistics.median(
                        q["elapsed_seconds"] for q in result["queries"])
                    reference = base / f"reports/cpp-{count}-query.tsv"
                    result["outputs_match_cpp"] = (sha(output) == sha(reference)) if reference.exists() else None
                measure_bulk(base, kind, count, args.threads, bulk, result)
            results.append(result)
            save(results_path, results)
            print(json.dumps(result), flush=True)
            if result["build"]["exit_code"] != 0:
                raise RuntimeError(f"Build failed for {kind}-{count}; larger cases are stopped")
        oracle = base / "bin/oracle-real"
        if oracle.exists() and not any(flag in args.index_flags for flag in ["--log-count", "--log-base", "--log-bits", "--nocount", "--quantization"]):
            oracle_result = measure([oracle, query, fof, base / "reports/oracle", args.threads],
                                    base / f"logs/oracle-{count}.log")
            save(base / f"reports/oracle-timing-{count}.json", oracle_result)
            if oracle_result["exit_code"] != 0:
                raise RuntimeError("Independent FASTA oracle failed")
            successful = [kind for kind in args.kinds if any(r["kind"] == kind and r["datasets"] == count
                          and r["build"]["exit_code"] == 0 and r.get("queries")
                          and all(q["exit_code"] == 0 for q in r["queries"]) for r in results)]
            validate(base, count, successful)
        count *= 2
    (base / ("reports/complete-" + "_".join(args.kinds))).touch()


if __name__ == "__main__":
    main()
