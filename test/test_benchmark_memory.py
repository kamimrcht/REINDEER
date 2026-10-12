#!/usr/bin/env python3
"""Check the benchmark guard without allocating RAM on the shared server."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from benchmark_server import memory_headroom, memory_emergency

GIB = 1024**3
with tempfile.TemporaryDirectory() as temporary:
    stats = Path(temporary) / "arcstats"

    def arc(size, minimum, mru, mfu, extra=""):
        stats.write_text("name type data\n" + "\n".join(
            f"{name} 4 {value * GIB}" for name, value in
            [("size", size), ("c_min", minimum),
             ("mru_evictable_data", mru), ("mfu_evictable_data", mfu)]) + "\n" + extra)

    assert memory_headroom(20 * GIB, stats)["effective_available_bytes"] == 20 * GIB
    arc(240, 16, 140, 70)
    headroom = memory_headroom(20 * GIB, stats)
    assert headroom["effective_available_bytes"] == 125 * GIB
    assert not memory_emergency(headroom)
    assert memory_emergency(memory_headroom(7 * GIB, stats))
    arc(40, 16, 140, 70)  # Cap credit above the ARC minimum, even on racing reads.
    assert memory_headroom(20 * GIB, stats)["zfs_credit_bytes"] == 12 * GIB
    arc(16, 16, 140, 70)
    assert memory_headroom(20 * GIB, stats)["zfs_credit_bytes"] == 0
    arc(240, 16, 0, 0, "mru_ghost_evictable_data 4 999999999999\nmetadata_size 4 999999999999\n")
    assert memory_headroom(20 * GIB, stats)["zfs_credit_bytes"] == 0
    for text in ["size 4 bad\n", "size 4 -1\nc_min 4 0\nmru_evictable_data 4 0\nmfu_evictable_data 4 0\n"]:
        stats.write_text(text)
        assert memory_headroom(20 * GIB, stats)["zfs_credit_bytes"] == 0
    arc(240, 16, 140, 70)
    assert not memory_emergency(memory_headroom(8 * GIB, stats))
print("Benchmark memory guard regressions passed.")
