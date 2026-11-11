#!/bin/bash
#SBATCH --job-name=a2_partA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=data/slurm_%x_%j.out
# Add your cluster's partition here if required, e.g.:
##SBATCH --partition=<partition-name>

# Usage (from the archive root):
#   sbatch scripts/run_partA.sh sweep          # phase 1: tile sweep, find b*
#   sbatch scripts/run_partA.sh main <bstar>   # phase 2: everything else

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p data

PHASE=${1:?"phase required: sweep | main"}
REPS=10
BIN=./build/matmul_pthreads

# Record environment for the report
{
    echo "== date =="; date
    echo "== host =="; hostname
    echo "== compiler =="; g++ --version | head -1
    echo "== lscpu =="; lscpu
} > data/partA_env.txt

make clean >/dev/null
make

run() {  # args: M K N pr pc b outfile
    echo ">> $*"
    $BIN "$1" "$2" "$3" "$4" "$5" "$6" "$REPS" | tee -a "$7"
}

if [ "$PHASE" = "sweep" ]; then
    OUT=data/partA_sweep.txt
    : > "$OUT"
    for b in 0 16 32 64 128 256; do
        run 2048 2048 2048 0 0 "$b" "$OUT"
    done
    run 1024 1024 1024 0 0 0 "$OUT"
    echo "Sweep done. Pick b* = smallest median at 2048, then:"
    echo "  sbatch scripts/run_partA.sh main <bstar>"

elif [ "$PHASE" = "main" ]; then
    BSTAR=${2:?"b* required: sbatch scripts/run_partA.sh main <bstar>"}

    # --- One-dimensional partition: pr = p, pc = 1 ---
    OUT=data/partA_1d.txt; : > "$OUT"
    for n in 1024 2048; do
        for b in 0 "$BSTAR"; do
            for p in 1 2 4 8 12 18 36; do
                run $n $n $n "$p" 1 "$b" "$OUT"
            done
        done
    done

    # --- Two-dimensional square grids, b = b* ---
    OUT=data/partA_2d.txt; : > "$OUT"
    GRIDS="1x2 2x2 2x4 3x4 3x6 6x6"
    for n in 1024 2048; do
        for g in $GRIDS; do
            run $n $n $n "${g%x*}" "${g#*x}" "$BSTAR" "$OUT"
        done
    done
    run 2048 2048 2048 1 36 "$BSTAR" "$OUT"     # column grid 1 x 36

    # --- Non-divisible dimensions ---
    OUT=data/partA_nondiv.txt; : > "$OUT"
    run 1000 1024 1536 6 6 "$BSTAR" "$OUT"

else
    echo "unknown phase: $PHASE" >&2
    exit 1
fi

# Collect RESULT lines into CSV for plotting
for f in data/partA_*.txt; do
    case "$f" in *env*) continue ;; esac
    grep '^RESULT' "$f" | sed 's/RESULT //; s/[A-Za-z_]*=//g; s/ /,/g' \
        > "${f%.txt}.csv" || true
done