#!/usr/bin/env python3
"""Split the 1,529 pairs into K record-shards for a parallel search array.

The serial design (all records in one DB, one searchmatrix) threw away the
per-chunk parallelism the real sweep relies on: ~14 s/record * 6,336 queries is
~6 h for one side alone. Sharding the records across array tasks restores it --
each task builds its own small exp+pred DB and searches both, so wall time is
(records_per_shard * 14 s * 2 sides) plus a ~1 min build.

Both halves of a given pair stay in the SAME shard so a task is self-contained;
which shard a pair lands in does not affect the per-domain result (each record is
searched independently), so a plain round-robin is fine.
"""
import json, sys
from pathlib import Path

R = Path("/home/rschaeff/work/prosmos_2026/instr_agree")
K = int(sys.argv[1]) if len(sys.argv) > 1 else 32
pairs = json.load(open(R / "pairs.json"))
shard_dir = R / "shards"
shard_dir.mkdir(exist_ok=True)
for old in shard_dir.glob("*.txt"):
    old.unlink()

for i in range(K):
    mine = pairs[i::K]
    with open(shard_dir / f"exp_{i:03d}.txt", "w") as fh:
        fh.write("".join(p["exp_pdb"] + "\n" for p in mine))
    with open(shard_dir / f"pred_{i:03d}.txt", "w") as fh:
        fh.write("".join(p["pred_pdb"] + "\n" for p in mine))
print(f"wrote {K} shards, ~{len(pairs)//K} pairs each, to {shard_dir}")
print(f"K={K}")
