#!/usr/bin/env python3
"""Write a live, auditable CSV/JSON/text analysis of the server benchmark."""
import argparse
import csv
import datetime
import io
import json
import os
from pathlib import Path
import re
import time


def replace(path, content):
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(content)
    temporary.replace(path)


def read(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default


def summarize(base, results_file):
    rows = read(results_file, [])
    metadata = read(base / "reports/environment.json", {})
    complete = base / ("reports/complete-" + results_file.stem.removeprefix("results-"))
    details = []
    for row in rows:
        count, kind, build = row["datasets"], row["kind"], row["build"]
        dedicated = base / f"reports/correctness-{count}-{kind}.json"
        correctness = read(dedicated, read(base / f"reports/correctness-{count}.json", {}))
        oracle = correctness.get("implementations", {}).get(kind, {})
        short = row.get("queries", [])
        bulk = row.get("bulk_queries", [])
        sample = {
            "datasets": count, "kind": kind, "threads": row["threads"],
            "binary_sha256": row["binary_sha256"], "build_exit_code": build["exit_code"],
            "input_plain_gb": row.get("input_plain_bytes", 0) / 10**9,
            "input_compressed_gb": row.get("input_compressed_bytes", 0) / 10**9,
            "build_seconds": build.get("elapsed_seconds"),
            "build_peak_rss_gib": build.get("peak_rss_kib", 0) / 1024**2,
            "build_cpu_cores_average": (build.get("user_seconds", 0) + build.get("system_seconds", 0)) / build["elapsed_seconds"]
                if build.get("elapsed_seconds", 0) else None,
            "index_gb": row.get("index_bytes", 0) / 10**9,
            "construction_disk_peak_gb_sampled": build.get("peak_construction_disk_bytes_sampled", 0) / 10**9,
            "short_query_seconds_median": row.get("query_seconds_median"),
            "bulk_query_seconds_median": row.get("bulk_query_seconds_median"),
            "query_peak_rss_gib": max((q.get("peak_rss_kib", 0) for q in short + bulk), default=0) / 1024**2,
            "oracle_query_dataset_pairs": oracle.get("query_dataset_pairs_checked"),
            "oracle_mismatches": oracle.get("mismatches"),
            "short_outputs_match_cpp": row.get("outputs_match_cpp"),
            "bulk_outputs_match_cpp": row.get("bulk_outputs_match_cpp"),
            "query_outputs_repeatable": row.get("outputs_repeatable") and row.get("bulk_outputs_repeatable"),
            "resource_limit": build.get("resource_limit"),
            "index_flags": build["command"][build["command"].index("-t") + 2:],
            "temporary_roots_measured": build.get("temporary_roots_measured", []),
            "index_components_bytes": row.get("index_components_bytes", {}),
        }
        log = base / f"logs/{kind}-{count}-build.log"
        if log.exists():
            text = log.read_text(errors="replace")
            match = re.search(r"SSHash construction: ([\d.]+)s; (.*)", text)
            if match:
                sample["sshash_seconds"] = float(match[1])
                sample["sshash_components_bytes"] = {key: int(value) for key, value in re.findall(r"(\w+)=(\d+)B", match[2])}
            match = re.search(r"Union graph: (\d+) abundance labels, (\d+) strings, (\d+) shared profiles", text)
            if match:
                sample.update(zip(["union_labels", "union_strings", "union_profiles"], map(int, match.groups())))
            sample["ggcat_phase_seconds"] = {name: float(seconds) for name, seconds in re.findall(
                r"Finished phase: ([^.]+)\. phase duration: ([\d.]+)s", text)}
        details.append(sample)
    query = read(base / "reports/queries-throughput.json", {})
    for sample in details:
        if sample["bulk_query_seconds_median"]:
            sample["bulk_kmers_per_second_including_load"] = query.get("kmers", 0) / sample["bulk_query_seconds_median"]
    comparisons = []
    for sample in details:
        reference = next((r for r in details if r["datasets"] == sample["datasets"] and r["kind"] == "cpp"), None)
        if sample["kind"] == "cpp" or not reference or reference["build_exit_code"] or sample["build_exit_code"]:
            continue
        comparison = {"datasets": sample["datasets"], "kind": sample["kind"]}
        for key, label in [("build_seconds", "build_speedup"), ("bulk_query_seconds_median", "bulk_query_speedup"),
                           ("index_gb", "index_size_ratio"), ("build_peak_rss_gib", "construction_ram_ratio"),
                           ("construction_disk_peak_gb_sampled", "construction_disk_ratio")]:
            a, b = reference.get(key), sample.get(key)
            if a and b:
                comparison[label] = a / b if "speedup" in label else b / a
        comparisons.append(comparison)
    run_config = read(base / ("reports/run-config-" + results_file.stem.removeprefix("results-") + ".json"), {})
    maximum = run_config.get("maximum", metadata.get("maximum_datasets", 1024))
    minimum = run_config.get("minimum", 2)
    expected_counts = []
    count = minimum
    while count <= maximum:
        expected_counts.append(count)
        count *= 2
    kinds = set(results_file.stem.removeprefix("results-").split("_"))
    verified_counts = [count for count in expected_counts if kinds and all(
        any(r["datasets"] == count and r["kind"] == kind and r["build_exit_code"] == 0
            and r["oracle_mismatches"] == 0 and r["query_outputs_repeatable"]
            and r["short_outputs_match_cpp"] is not False
            and r["bulk_outputs_match_cpp"] is not False for r in details)
        for kind in kinds)]
    verified_comparisons = [count for count in verified_counts if
        any(r["datasets"] == count and r["kind"] == "cpp" and r["build_exit_code"] == 0
            and r["oracle_mismatches"] == 0 for r in details)
        and any(r["datasets"] == count and r["kind"] != "cpp" and r["build_exit_code"] == 0
            and r["oracle_mismatches"] == 0 and r["short_outputs_match_cpp"]
            and r["bulk_outputs_match_cpp"] for r in details)]
    result = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "benchmark_directory": str(base), "results_file": str(results_file),
        "completed": complete.exists(), "all_cases_attempted": complete.exists(),
        "all_requested_cases_succeeded": complete.exists() and len(verified_counts) == len(expected_counts),
        "requested_maximum": maximum,
        "requested_minimum": minimum,
        "largest_independently_verified_case": max(verified_counts, default=0),
        "largest_independently_verified_comparison": max(verified_comparisons, default=0),
        "largest_case_recorded": max((r["datasets"] for r in details), default=0),
        "environment": metadata, "query_workload": query,
        "method": "Nested deterministic dataset prefixes; k=31; 64 threads; native CPUs; query median of three fresh processes, including loading, with warm OS cache; source abundance oracle for each prefix. Implementations build from identical plain unitigs wherever compared; staging/decompression is recorded separately and excluded from construction times. Requested implementations: " + ", ".join(sorted(kinds)) + ".",
        "limitations": [
            "Server is shared; no global cache flush or NUMA isolation.",
            "GNU time peak RSS covers process memory; it does not attribute filesystem cache to the process.",
            "Construction disk is sampled every 20 seconds and is a lower bound on peak; early retained C++ cases have legacy measurement roots.",
            "C++ includes file-descriptor and nested-OpenMP runtime corrections, with the original indexing algorithm.",
            "Rust tuning is explicit in index_flags; legacy CLI defaults are unchanged.",
            "Source oracle covers the fixed query set, not every indexed k-mer; bulk outputs are also compared byte-for-byte wherever a completed C++ reference exists.",
            "Missing/failed cases are not performance results. A run is complete only when the completion marker exists.",
        ],
        "measurements": details, "comparisons": comparisons,
    }
    return result


def write(base, results_file):
    result = summarize(base, results_file)
    prefix = base / ("reports/analysis-" + results_file.stem.removeprefix("results-"))
    replace(prefix.with_suffix(".json"), json.dumps(result, indent=2) + "\n")
    fields = ["datasets", "kind", "threads", "build_exit_code", "build_seconds", "build_peak_rss_gib",
              "build_cpu_cores_average", "index_gb", "construction_disk_peak_gb_sampled", "short_query_seconds_median",
              "bulk_query_seconds_median", "bulk_kmers_per_second_including_load", "query_peak_rss_gib",
              "oracle_query_dataset_pairs", "oracle_mismatches", "short_outputs_match_cpp", "bulk_outputs_match_cpp", "resource_limit"]
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(result["measurements"])
    replace(prefix.with_suffix(".csv"), stream.getvalue())
    lines = [f"Updated: {result['generated_utc']}",
             f"All cases attempted: {result['all_cases_attempted']}; all succeeded: {result['all_requested_cases_succeeded']}; requested maximum: {result['requested_maximum']}",
             f"Largest independently verified comparison: {result['largest_independently_verified_comparison']}",
             result["method"], "", "N  implementation  build_s  build_GiB  index_GB  peak_disk_GB  bulk_query_s  oracle_mismatches"]
    for r in result["measurements"]:
        lines.append(f"{r['datasets']:4} {r['kind']:16} {r['build_seconds']} {r['build_peak_rss_gib']:.3f} {r['index_gb']:.3f} "
                     f"{r['construction_disk_peak_gb_sampled']:.3f} {r['bulk_query_seconds_median']} {r['oracle_mismatches']}")
    lines += ["", "Ratios against C++:"] + [json.dumps(c) for c in result["comparisons"]]
    lines += ["", "Limitations:"] + result["limitations"]
    replace(prefix.with_suffix(".txt"), "\n".join(lines) + "\n")
    return result["completed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--results", default="results-cpp_rust-wide.json")
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    base = args.base.resolve()
    while True:
        complete = write(base, base / "reports" / args.results)
        if complete or not args.watch:
            break
        time.sleep(30)


if __name__ == "__main__":
    main()
