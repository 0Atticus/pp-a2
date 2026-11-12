#!/bin/bash
#SBATCH --job-name=a2_all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=slurm_all_%j.out
#SBATCH --partition=Centaurus

# Master script: ONE job, ONE node, so T_seq and every T(p) of Parts A and B are
# measured under identical conditions (the spec requires this).
#
# Order: build -> Part A tile sweep -> pick b* -> Part A experiments
#        -> Part B experiments -> collect everything into data/all_results.csv
#
# Usage (from the archive root, next to the Makefile):
#   sbatch scripts/run_all.sbatch                              # b* chosen automatically
#   sbatch --export=ALL,BSTAR=64 scripts/run_all.sbatch        # force a b*
#
# Dry run on a small problem (any machine, no Slurm):
#   BIG=256 SMALLN=128 NDM=100 NDK=128 NDN=150 REPS=2 bash scripts/run_all.sbatch

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p data

REPS=${REPS:-10}
BIG=${BIG:-2048}          # large matrix size
SMALLN=${SMALLN:-1024}    # small matrix size
NDM=${NDM:-1000}; NDK=${NDK:-1024}; NDN=${NDN:-1536}   # non-divisible case
BSTAR=${BSTAR:-}          # optional override of the automatic choice

PTH=./build/matmul_pthreads
OMP=./build/matmul_omp
GRIDS="1x2 2x2 2x4 3x4 3x6 6x6"
PLIST="1 2 4 8 12 18 36"

# ---------------------------------------------------------------- environment
{
    echo "== date =="; date
    echo "== host =="; hostname
    echo "== slurm job =="; echo "${SLURM_JOB_ID:-none}"
    echo "== compiler =="; g++ --version | head -1
    echo "== OMP env (Part B runs) =="; echo "OMP_PROC_BIND=true OMP_PLACES=cores"
    echo "== settings =="; echo "REPS=$REPS BIG=$BIG SMALLN=$SMALLN"
    echo "== lscpu =="; lscpu
} > data/env.txt
cp data/env.txt data/partA_env.txt
cp data/env.txt data/partB_env.txt

make clean >/dev/null
make all || { echo "build failed" >&2; exit 1; }

# -------------------------------------------------------------------- runners
runA() {  # M K N pr pc b outfile
    echo ">> A: $*"
    $PTH "$1" "$2" "$3" "$4" "$5" "$6" "$REPS" | tee -a "$7" \
        || echo "FAILED: A $*" | tee -a "$7"
}
runB() {  # M K N pr pc b outfile
    echo ">> B: $*"
    OMP_PROC_BIND=true OMP_PLACES=cores \
        $OMP "$1" "$2" "$3" "$4" "$5" "$6" "$REPS" | tee -a "$7" \
        || echo "FAILED: B $*" | tee -a "$7"
}

# ----------------------------------------------------------- Part A: tile sweep
OUT=data/partA_sweep.txt; : > "$OUT"
for b in 0 16 32 64 128 256; do
    runA $BIG $BIG $BIG 0 0 "$b" "$OUT"
done
runA $SMALLN $SMALLN $SMALLN 0 0 0 "$OUT"

# Pick b*: smallest median at the large size among verified runs
grep '^RESULT' "$OUT" | sed 's/RESULT //; s/[A-Za-z_]*=//g; s/ /,/g' > data/partA_sweep.csv
echo "== sweep medians at $BIG (sorted) =="
awk -F, -v n="$BIG" '$1==n && $14=="yes" {print $10, $7}' data/partA_sweep.csv | sort -g \
    | tee data/sweep_ranking.txt
AUTO=$(head -1 data/sweep_ranking.txt | awk '{print $2}')
if awk 'NR==1{a=$1} NR==2{b=$1} END{exit !(NR>=2 && (b-a)/a < 0.01)}' data/sweep_ranking.txt; then
    echo "WARNING: top two tile sizes are within 1% of each other; b* is essentially a tie."
fi
if [ -z "$BSTAR" ]; then BSTAR=$AUTO; fi
if [ -z "$BSTAR" ]; then
    echo "could not determine b* (no verified sweep runs)" >&2
    exit 1
fi
echo "b* = $BSTAR (automatic choice was $AUTO)" | tee data/bstar.txt

# ------------------------------------------------------------- Part A: main
OUT=data/partA_1d.txt; : > "$OUT"
for n in $SMALLN $BIG; do
    for b in 0 "$BSTAR"; do
        for p in $PLIST; do runA $n $n $n "$p" 1 "$b" "$OUT"; done
    done
done

OUT=data/partA_2d.txt; : > "$OUT"
for n in $SMALLN $BIG; do
    for g in $GRIDS; do runA $n $n $n "${g%x*}" "${g#*x}" "$BSTAR" "$OUT"; done
done
runA $BIG $BIG $BIG 1 36 "$BSTAR" "$OUT"          # column grid 1 x 36

OUT=data/partA_nondiv.txt; : > "$OUT"
runA $NDM $NDK $NDN 6 6 "$BSTAR" "$OUT"

# ------------------------------------------------------------- Part B
OUT=data/partB_seq.txt; : > "$OUT"
for n in $SMALLN $BIG; do
    for b in 0 "$BSTAR"; do runB $n $n $n 0 0 "$b" "$OUT"; done
done

OUT=data/partB_1d.txt; : > "$OUT"
for n in $SMALLN $BIG; do
    for b in 0 "$BSTAR"; do
        for p in $PLIST; do runB $n $n $n "$p" 1 "$b" "$OUT"; done
    done
done

OUT=data/partB_2d.txt; : > "$OUT"
for n in $SMALLN $BIG; do
    for g in $GRIDS; do runB $n $n $n "${g%x*}" "${g#*x}" "$BSTAR" "$OUT"; done
done
runB $BIG $BIG $BIG 1 36 "$BSTAR" "$OUT"

OUT=data/partB_nondiv.txt; : > "$OUT"
runB $NDM $NDK $NDN 6 6 "$BSTAR" "$OUT"

# ----------------------------------------------------------------- collect
bash scripts/collect.sh data

echo "== run counts (expected in parentheses) =="
chk() { printf "%-22s %3s (%s)\n" "$1" "$(grep -c '^RESULT' "data/$1.txt")" "$2"; }
chk partA_sweep 7; chk partA_1d 28; chk partA_2d 13; chk partA_nondiv 1
chk partB_seq 4;   chk partB_1d 28; chk partB_2d 13; chk partB_nondiv 1

echo "== FAILED or unverified runs =="
grep -h '^FAILED' data/part[AB]_*.txt || echo "  no FAILED lines"
grep -h 'verified=NO\|mapped=NO' data/part[AB]_*.txt || echo "  no verification failures"

tar czf "data_partAB_${SLURM_JOB_ID:-local}.tgz" data
echo "archived: data_partAB_${SLURM_JOB_ID:-local}.tgz"
