#!/usr/bin/env python3
"""Exercise the public CLI against an independently computed k-mer oracle."""
import csv
import math
import os
import random
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BINARY = Path(os.environ.get("REINDEER_BINARY", ROOT / "Reindeer")).resolve()
K = 31


def canonical(sequence):
    return min(sequence, sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])


def run(*args):
    result = subprocess.run([str(BINARY), *map(str, args)],
                            cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr


def raw(values):
    runs = []
    for position, value in enumerate(values):
        value = str(value) if value else "*"
        if runs and runs[-1][2] == value:
            runs[-1][1] = position
        else:
            runs.append([position, position, value])
    return ",".join(f"{start}-{end}:{value}" for start, end, value in runs)


def main():
    rng = random.Random(762)
    sequences = ["".join(rng.choices("ACGT", k=K)) for _ in range(15)]
    long_sequence = "".join(rng.choices("ACGT", k=K + 12))
    abundances = [1, 2, 3, 4, 27, 101, 40812, 65535, 65536, 0, 1.5, 2.5]
    queries = {f"count_{i}": sequences[i] for i in range(len(abundances))}
    queries["missing"] = sequences[-1]
    queries["partial_half"] = sequences[4] + "A"
    queries["partial_quarter"] = sequences[4] + "AAA"
    queries["repeated_monotig"] = long_sequence
    queries["reverse"] = long_sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    with tempfile.TemporaryDirectory(prefix="reindeer-modes-") as tmp:
        tmp = Path(tmp)
        samples = [tmp / f"sample{i}.fa" for i in range(3)]
        data = [list(zip(sequences, abundances)) + [(long_sequence, 27)],
                [(sequences[4], 2), (sequences[5], 4), (long_sequence, 8)],
                [(sequences[-2], 3)]]
        for sample, entries in zip(samples, data):
            # Include all four supported tags, at end of header and with links.
            sample.write_text("".join(
                f">{i} {['km', 'KM', 'ka', 'KA'][i % 4]}:f:{count}"
                f"{' L:+:0:+' if i % 2 else ''}\n{sequence}\n"
                for i, (sequence, count) in enumerate(entries)))
        fof = tmp / "fof.txt"
        fof.write_text("".join(f"{sample}\n" for sample in samples))
        query_file = tmp / "queries.fa"
        query_file.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in queries.items()))

        checks = 0
        for mode, flags in [("raw", []), ("log", ["--log-count"]),
                            ("presence", ["--nocount"])]:
            oracle = []
            for entries in data:
                counts = {}
                for sequence, count in entries:
                    value = (1 if mode == "presence" else
                             (max(1, math.ceil(math.log2(count))) if count else 0)
                             if mode == "log" else int(min(count, 65535)))
                    for i in range(len(sequence) - K + 1):
                        counts[canonical(sequence[i:i + K])] = value
                oracle.append(counts)
            for memory in [False, True]:
                index = tmp / f"{mode}-{memory}"
                run("--index", "-f", fof, "-o", index, "-t", 2, *flags,
                    *(["--mem-query"] if memory else []))
                for threshold in [0, 40, 50, 51, 100]:
                    for fmt in ["raw", "sum", "average", "mean"]:
                        output = tmp / "result.tsv"
                        # Mode comes from metadata; no query-side mode flag.
                        run("--query", "-l", index, "-q", query_file,
                            "-o", output, "-P", threshold, "--format", fmt,
                            "-t", 2)
                        with output.open() as stream:
                            rows = list(csv.reader(stream, delimiter="\t"))
                        assert rows[0] == ["query", "sample0", "sample1", "sample2"], rows[0]
                        assert [row[0] for row in rows[1:]] == list(queries)
                        for row in rows[1:]:
                            sequence = queries[row[0]]
                            for sample, actual in enumerate(row[1:]):
                                values = [oracle[sample].get(canonical(sequence[i:i + K]), 0)
                                          for i in range(len(sequence) - K + 1)]
                                coverage = sum(value > 0 for value in values) * 100 // len(values)
                                if coverage < threshold:
                                    expected = "*"
                                elif fmt == "raw":
                                    expected = raw(values)
                                elif fmt == "sum":
                                    expected = str(sum(values))
                                else:
                                    expected = f"{sum(values) / len(values):.2f}"
                                assert actual == expected, (mode, memory, threshold, fmt, row[0], sample, actual, expected)
                                checks += 1
        print(f"Mode regressions passed ({checks} assertions, disk and memory indexes).")

        bcalm = ROOT / "bin/bcalm"
        if bcalm.exists():
            reads = ROOT / "test/input_bcalm/SRR10092187_10k.fastq"
            prefix = tmp / "bcalm-smoke"
            result = subprocess.run(
                [str(bcalm), "-in", str(reads), "-kmer-size", str(K),
                 "-abundance-min", "2", "-nb-cores", "2", "-out", str(prefix)],
                cwd=tmp, capture_output=True, text=True, timeout=120)
            assert result.returncode == 0, result.stdout + result.stderr
            frequencies = Counter(
                canonical(sequence[i:i + K])
                for sequence in reads.read_text().splitlines()[1::4]
                for i in range(len(sequence) - K + 1)
                if set(sequence[i:i + K]) <= set("ACGT"))
            expected = {sequence for sequence, count in frequencies.items() if count >= 2}
            actual = {canonical(sequence[i:i + K])
                      for sequence in Path(str(prefix) + ".unitigs.fa").read_text().splitlines()[1::2]
                      for i in range(len(sequence) - K + 1)}
            assert actual == expected, "BCALM changed the solid k-mer spectrum"
            print(f"BCALM regression passed ({len(expected)} solid k-mers checked against reads).")


if __name__ == "__main__":
    main()
