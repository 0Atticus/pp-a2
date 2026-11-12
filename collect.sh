#!/bin/bash
# Merge every RESULT line in data/part[AB]_*.txt into data/all_results.csv and
# compute S = T_seq / median, where T_seq is the median of the untiled (b=0)
# sequential run of the SAME part and SAME size. Safe to re-run any time.
#
# Usage: bash scripts/collect.sh [datadir]     (default: data)

set -uo pipefail
D=${1:-data}
TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT

for f in "$D"/partA_*.txt "$D"/partB_*.txt; do
    [ -e "$f" ] || continue
    base=$(basename "$f" .txt)
    case "$base" in *env*) continue ;; esac
    part=${base%%_*}        # partA | partB
    exp=${base#*_}          # sweep | seq | 1d | 2d | nondiv
    grep '^RESULT' "$f" | sed "s/^RESULT /$part $exp /" >> "$TMP"
done

if [ ! -s "$TMP" ]; then
    echo "no RESULT lines found in $D" >&2
    exit 1
fi

awk '
function parse(   i, kv) {
    delete v
    for (i = 3; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
}
# pass 1: collect T_seq (untiled sequential medians) per part and size
FNR == NR {
    parse()
    if (v["pr"] == 0 && v["pc"] == 0 && v["b"] == 0)
        tseq[$1 SUBSEP v["M"] SUBSEP v["K"] SUBSEP v["N"]] = v["median"]
    next
}
# pass 2: print rows
FNR == 1 {
    print "part,experiment,M,K,N,pr,pc,p,b,reps,avg,median,min,max,gflops,speedup,verified,mapped"
}
{
    parse()
    key = $1 SUBSEP v["M"] SUBSEP v["K"] SUBSEP v["N"]
    sp = ""
    if ((key in tseq) && v["median"] > 0) sp = sprintf("%.4f", tseq[key] / v["median"])
    mp = ("mapped" in v) ? v["mapped"] : "na"
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        $1, $2, v["M"], v["K"], v["N"], v["pr"], v["pc"], v["p"], v["b"], v["reps"],
        v["avg"], v["median"], v["min"], v["max"], v["gflops"], sp, v["verified"], mp
}
' "$TMP" "$TMP" > "$D/all_results.csv"

rows=$(($(wc -l < "$D/all_results.csv") - 1))
bad=$(grep -c -E ',NO(,|$)' "$D/all_results.csv" || true)
echo "wrote $D/all_results.csv: $rows rows, $bad with verified=NO or mapped=NO"
