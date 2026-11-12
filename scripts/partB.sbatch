#!/bin/bash
#SBATCH --job-name=a2_partB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=slurm_partB_%j.out
#SBATCH --partition=Centaurus

# Usage (from the archive root, next to the Makefile):
#   sbatch scripts/run_partB.sbatch <bstar>      e.g. sbatch scripts/run_partB.sbatch 64
# Use the same b* as Part A.

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p data

BSTAR=${1:?"b* required: sbatch scripts/run_partB.sbatch <bstar>"}
REPS=10
BIN=./build/matmul_omp

# Required by the assignment for Part B
export OMP_PROC_BIND=true
export OMP_PLACES=cores

{
    echo "== date =="; date
    echo "== host =="; hostname
    echo "== compiler =="; g++ --version | head -1
    echo "== OMP env =="; echo "OMP_PROC_BIND=$OMP_PROC_BIND OMP_PLACES=$OMP_PLACES"
    echo "== lscpu =="; lscpu
} > data/partB_env.txt

make clean >/dev/null
make omp || { echo "build failed" >&2; exit 1; }

run() {  # args: M K N pr pc b outfile
    echo ">> $*"
    $BIN "$1" "$2" "$3" "$4" "$5" "$6" "$REPS" | tee -a "$7" \
        || echo "FAILED: $*" | tee -a "$7"
}

# --- Sequential baselines (b = 0 and b = b*) at both sizes ---
OUT=data/partB_seq.txt; : > "$OUT"
for n in 1024 2048; do
    for b in 0 "$BSTAR"; do
        run $n $n $n 0 0 "$b" "$OUT"
    done
done

# --- One-dimensional partition: pr = p, pc = 1 ---
OUT=data/partB_1d.txt; : > "$OUT"
for n in 1024 2048; do
    for b in 0 "$BSTAR"; do
        for p in 1 2 4 8 12 18 36; do
            run $n $n $n "$p" 1 "$b" "$OUT"
        done
    done
done

# --- Two-dimensional square grids, b = b* ---
OUT=data/partB_2d.txt; : > "$OUT"
GRIDS="1x2 2x2 2x4 3x4 3x6 6x6"
for n in 1024 2048; do
    for g in $GRIDS; do
        run $n $n $n "${g%x*}" "${g#*x}" "$BSTAR" "$OUT"
    done
done
run 2048 2048 2048 1 36 "$BSTAR" "$OUT"     # column grid 1 x 36

# --- Non-divisible dimensions ---
OUT=data/partB_nondiv.txt; : > "$OUT"
run 1000 1024 1536 6 6 "$BSTAR" "$OUT"

# Collect RESULT lines into CSV (columns: M,K,N,pr,pc,p,b,reps,avg,median,min,max,gflops,verified,mapped)
for f in data/partB_*.txt; do
    case "$f" in *env*) continue ;; esac
    grep '^RESULT' "$f" | sed 's/RESULT //; s/[A-Za-z_]*=//g; s/ /,/g' \
        > "${f%.txt}.csv" || true
done

echo "FAILED lines:"; grep -h '^FAILED' data/partB_*.txt || echo "  none"