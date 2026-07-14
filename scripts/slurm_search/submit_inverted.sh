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
# local: each task writes searchmatrix's one-file-per-hit output to the compute
# node's own disk and ships back ONE compressed TSV. tree: the old layout, tiny
# files straight onto NFS (~24M inodes for a full AFDB sweep). Default local.
: "${HITS_MODE:=local}"
[ -d "$QDIR" ] || { echo "QDIR not a dir: $QDIR" >&2; exit 1; }
[ -f "$DB" ]   || { echo "DB not found: $DB" >&2; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$OUT/logs" "$OUT/parts" "$OUT/chunks"
echo "$HITS_MODE" > "$OUT/.hits_mode"   # retry_missing.sh must match the run's layout

# 1) query manifest (all queries, one absolute path per line)
MANIFEST="$OUT/manifest.txt"
find "$QDIR" -name '*.query' -type f | sort > "$MANIFEST"
NQ=$(wc -l < "$MANIFEST"); [ "$NQ" -gt 0 ] || { echo "no queries in $QDIR" >&2; exit 1; }
echo "queries: $NQ  hits mode: $HITS_MODE" >&2

# 2) tree mode only: pre-create the per-query NFS dirs (searchmatrix writes into
# them but does not create them). In local mode each task makes them on its own
# node instead, so nothing but the per-chunk TSV ever lands on NFS.
if [ "$HITS_MODE" = "tree" ]; then
    mkdir -p "$OUT/hits"
    echo "pre-creating $NQ hit dirs..." >&2
    while read -r q; do n=$(basename "$q" .query); mkdir -p "$OUT/hits/$n"; done < "$MANIFEST"
else
    mkdir -p "$OUT/hitparts"
fi

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
        --export=CHUNK_LIST="$CHUNK_LIST",MANIFEST="$MANIFEST",HITS_ROOT="$OUT/hits/",OUT="$OUT",OFFSET="$offset",HITS_MODE="$HITS_MODE",ARCHIVE_PY="$HERE/archive_hits.py" \
        "$HERE/array_inverted.sbatch")
    ARRAY_IDS+=("$JID"); echo "  array $JID (offset $offset, $this tasks${PREV:+, after $PREV})" >&2
    PREV="$JID"; offset=$((offset + this))
done
DEPS=$(IFS=:; echo "${ARRAY_IDS[*]}")
# Summary only. Deliberately NOT `find $OUT/hits -name '*.txt'` -- the tree holds
# tens of millions of files and walking it over NFS costs hours, and `sort -u` on
# the basenames discards which QUERY hit which record. submit_archive.sh folds the
# tree into one zstd table and emits distinct_hitters.txt from that single pass.
MID=$(sbatch --parsable --job-name=prosmos-inv-merge --time=00:20:00 --mem=1G \
    --dependency=afterany:${DEPS} --chdir="$OUT/logs" \
    --wrap="{ printf 'chunk\truntime_sec\texit_code\trecords\n'; cat $OUT/parts/*.tsv 2>/dev/null | sort; } > $OUT/summary.tsv; echo done")

# In local mode the archive is just a concatenation of the per-chunk TSVs -- no
# tree walk, no separate submit_archive.sh pass. merge_chunks refuses to build it
# if any chunk is missing a TSV or left a .partial behind.
if [ "$HITS_MODE" = "local" ]; then
    PY=/sw/apps/Anaconda3-2023.09-0/bin/python
    AID=$(sbatch --parsable --job-name=prosmos-inv-archive --time=02:00:00 --mem=8G \
        --dependency=afterok:${MID} --chdir="$OUT/logs" \
        --wrap="$PY $HERE/archive_hits.py merge_chunks $OUT")
    echo "archive job: $AID" >&2
fi
echo "merge job: $MID" >&2
echo "watch: squeue -j ${ARRAY_IDS[*]}" >&2
