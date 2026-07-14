#!/bin/bash
# Resubmit the chunks of a loop-inverted sweep that never produced a clean part
# file (TIMEOUT, node failure, rc!=0). Chunk runtimes have a heavy tail -- on the
# PDB DB the slowest chunks are 3-4x the mean -- so the original 8h --time kills a
# handful of them. SLURM does not let a non-admin RAISE --time on an already
# submitted job, so the only fix is to resubmit the stragglers with a bigger limit.
#
# Reruns are safe/idempotent: searchmatrix writes one hit file per DB record
# (HITS_ROOT/<query>/<recordid>.txt), so redoing a chunk rewrites the same paths.
#
# Usage:
#   OUT=path/to/run [TIMELIMIT=24:00:00] [CONCURRENCY=32] ./retry_missing.sh
set -euo pipefail
: "${OUT:?missing OUT}"
: "${TIMELIMIT:=24:00:00}"
: "${CONCURRENCY:=32}"
HERE=$(cd "$(dirname "$0")" && pwd)

CHUNK_LIST="$OUT/chunks.list"
MANIFEST="$OUT/manifest.txt"
[ -f "$CHUNK_LIST" ] || { echo "no $CHUNK_LIST" >&2; exit 1; }
[ -f "$MANIFEST" ]   || { echo "no $MANIFEST" >&2; exit 1; }

N=$(wc -l < "$CHUNK_LIST")

# A chunk is "done" iff SOME part file records its basename with exit_code 0.
# Keying on the chunk NAME (not the part's line number) keeps this correct across
# repeated retry passes, whose parts are numbered 1..n_missing, not 1..N.
RETRY_LIST="$OUT/chunks.retry.list"
python3 - "$CHUNK_LIST" "$OUT/parts" "$RETRY_LIST" <<'PY'
import sys, glob, os
chunk_list, parts_dir, out = sys.argv[1:4]
done = set()
for p in glob.glob(os.path.join(parts_dir, "*.tsv")):
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[2] == "0":
            done.add(f[0])
chunks = [l.strip() for l in open(chunk_list) if l.strip()]
missing = [c for c in chunks if os.path.basename(c) not in done]
open(out, "w").write("".join(c + "\n" for c in missing))
print(f"{len(missing)} missing of {len(chunks)}", file=sys.stderr)
PY

NM=$(wc -l < "$RETRY_LIST")
echo "missing/failed chunks: $NM / $N" >&2
[ "$NM" -gt 0 ] || { echo "nothing to retry" >&2; exit 0; }
JID=$(sbatch --parsable --time="$TIMELIMIT" \
    --array=1-${NM}%${CONCURRENCY} --chdir="$OUT/logs" \
    --job-name=prosmos-retry \
    --export=CHUNK_LIST="$RETRY_LIST",MANIFEST="$MANIFEST",HITS_ROOT="$OUT/hits/",OUT="$OUT",OFFSET=0,PARTS_PREFIX=retry \
    "$HERE/array_inverted.sbatch")
echo "retry array: $JID ($NM tasks, --time=$TIMELIMIT)" >&2

# Summary only. Deliberately NOT `find $OUT/hits -name '*.txt'` -- the tree holds
# tens of millions of files and walking it over NFS costs hours. submit_archive.sh
# emits distinct_hitters.txt (and the query->record mapping the find threw away)
# from the single pass it has to make anyway.
MID=$(sbatch --parsable --job-name=prosmos-retry-merge --time=00:20:00 --mem=1G \
    --dependency=afterany:"$JID" --chdir="$OUT/logs" \
    --wrap="{ printf 'chunk\truntime_sec\texit_code\trecords\n'; cat $OUT/parts/*.tsv 2>/dev/null | sort -u; } > $OUT/summary.tsv; echo done")
echo "retry merge: $MID" >&2
echo "when this and the sweep are done, archive the hits tree:" >&2
echo "  OUT=$OUT $HERE/submit_archive.sh" >&2
