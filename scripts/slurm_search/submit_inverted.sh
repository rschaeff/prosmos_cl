#!/bin/bash
# Submit a loop-inverted ProSMoS sweep: split the DB into chunks and run ALL
# queries (a manifest) against each chunk in one parse. Parses the DB once total
# instead of once per query (~2x less compute than submit.sh's query-per-task).
#
# Usage:
#   QDIR=path/to/queries DB=path/to/DB OUT=path/to/results \
#   [NCHUNK=800] [CONCURRENCY=200] [TIMELIMIT=08:00:00] ./submit_inverted.sh
set -euo pipefail
: "${QDIR:?missing QDIR}"; : "${DB:?missing DB}"; : "${OUT:?missing OUT}"
: "${NCHUNK:=800}"; : "${CONCURRENCY:=200}"; : "${TIMELIMIT:=08:00:00}"
[ -d "$QDIR" ] || { echo "QDIR not a dir: $QDIR" >&2; exit 1; }
[ -f "$DB" ]   || { echo "DB not found: $DB" >&2; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$OUT/logs" "$OUT/parts" "$OUT/chunks" "$OUT/hits"

# 1) query manifest (all queries, one absolute path per line)
MANIFEST="$OUT/manifest.txt"
find "$QDIR" -name '*.query' -type f | sort > "$MANIFEST"
NQ=$(wc -l < "$MANIFEST"); [ "$NQ" -gt 0 ] || { echo "no queries in $QDIR" >&2; exit 1; }
echo "queries: $NQ" >&2

# 2) pre-create per-query output dirs once (avoids each task's mkdir storm on NFS)
echo "pre-creating $NQ hit dirs..." >&2
while read -r q; do n=$(basename "$q" .query); mkdir -p "$OUT/hits/$n"; done < "$MANIFEST"

# 3) split the DB into NCHUNK chunks on record boundaries
echo "splitting DB into $NCHUNK chunks..." >&2
python3 - "$DB" "$OUT/chunks" "$NCHUNK" <<'PY'
import sys, re
db, outdir, nchunk = sys.argv[1], sys.argv[2], int(sys.argv[3])
NAME = re.compile(r'^\S+\.ssd\s')
total = sum(1 for line in open(db) if NAME.match(line))
per = total // nchunk + 1
ci = 0; rec = 0
f = open(f"{outdir}/chunk_{ci:04d}.db", "w"); cur = []
with open(db) as src:
    for line in src:
        if NAME.match(line):
            rec += 1
            if rec > 1 and rec % per == 1:
                f.write(''.join(cur)); cur = []; f.close()
                ci += 1; f = open(f"{outdir}/chunk_{ci:04d}.db", "w")
        cur.append(line)
f.write(''.join(cur)); f.close()
print(f"{total} records -> {ci+1} chunks (~{per}/chunk)")
PY
CHUNK_LIST="$OUT/chunks.list"
ls "$OUT/chunks"/chunk_*.db | sort > "$CHUNK_LIST"
N=$(wc -l < "$CHUNK_LIST")
echo "chunks: $N" >&2

# 4) submit array (chained chunks of <=MaxArraySize-1 to cap total concurrency)
MAX_ARRAY=$(scontrol show config 2>/dev/null | awk '/MaxArraySize/{print $3}'); : "${MAX_ARRAY:=1000}"
CH=$((MAX_ARRAY - 1))
ARRAY_IDS=(); PREV=""; offset=0
while [ "$offset" -lt "$N" ]; do
    remaining=$((N - offset)); this=$((remaining < CH ? remaining : CH))
    DEP=(); [ -n "$PREV" ] && DEP+=(--dependency=afterany:"$PREV")
    JID=$(sbatch --parsable --time="$TIMELIMIT" "${DEP[@]}" \
        --array=1-${this}%${CONCURRENCY} --chdir="$OUT/logs" \
        --export=CHUNK_LIST="$CHUNK_LIST",MANIFEST="$MANIFEST",HITS_ROOT="$OUT/hits/",OUT="$OUT",OFFSET="$offset" \
        "$HERE/array_inverted.sbatch")
    ARRAY_IDS+=("$JID"); echo "  array $JID (offset $offset, $this tasks${PREV:+, after $PREV})" >&2
    PREV="$JID"; offset=$((offset + this))
done
DEPS=$(IFS=:; echo "${ARRAY_IDS[*]}")
MID=$(sbatch --parsable --job-name=prosmos-inv-merge --time=00:20:00 --mem=1G \
    --dependency=afterany:${DEPS} --chdir="$OUT/logs" \
    --wrap="{ printf 'chunk\truntime_sec\texit_code\trecords\n'; cat $OUT/parts/*.tsv 2>/dev/null | sort; } > $OUT/summary.tsv; find $OUT/hits -name '*.txt' -printf '%f\n' | sed 's/\.txt\$//' | sort -u > $OUT/distinct_hitters.txt; echo done")
echo "merge job: $MID" >&2
echo "watch: squeue -j ${ARRAY_IDS[*]}" >&2
