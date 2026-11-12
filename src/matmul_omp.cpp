// ITCS 4145 Assignment 2, Part B: blocked matrix multiplication with OpenMP.
//
// Usage: ./matmul_omp M K N pr pc b [reps]
//   pr = pc = 0 : sequential implementation (no parallel region)
//   pr x pc     : block grid; ONE worksharing loop over the pr*pc blocks of C,
//                 schedule(static), run by p = pr*pc threads (iteration t -> thread t)
//   b           : cache tile size, 0 = no tiling
//   reps        : timed repetitions (default 10), after one untimed warm-up
//
// Constructs used: parallel for, schedule(static), shared clause,
// omp_get_thread_num, omp_get_num_threads, omp_set_num_threads, omp_get_wtime.
//
// Build: g++ -O2 -fopenmp matmul_omp.cpp -o matmul_omp
// Run with: OMP_PROC_BIND=true OMP_PLACES=cores

#include "assignment2.h"  // first, as the header asks

#include <omp.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// Copy a freshly generated matrix into caller-owned storage.
static void gen(double *X, int rows, int cols, uint64_t seed) {
    double *t = generate_matrix(rows, cols, seed);
    std::memcpy(X, t, (size_t)rows * cols * sizeof(double));
    free(t);
}

static inline int imin(int a, int b) { return a < b ? a : b; }

// Tiled kernel on block C[r0,r1) x [c0,c1). b == 0 means one tile (no tiling).
// Loop order inside a tile: i, k, j. Identical to Part A.
static void kernel(const double *A, const double *B, double *C,
                   int K, int N,
                   int r0, int r1, int c0, int c1, int b) {
    if (r1 <= r0 || c1 <= c0) return;
    const int tr = (b == 0) ? (r1 - r0) : b;
    const int tc = (b == 0) ? (c1 - c0) : b;
    const int tk = (b == 0) ? K : b;

    for (int ii = r0; ii < r1; ii += tr) {
        const int ie = imin(ii + tr, r1);
        for (int jj = c0; jj < c1; jj += tc) {
            const int je = imin(jj + tc, c1);

            for (int i = ii; i < ie; i++)
                for (int j = jj; j < je; j++) C[(size_t)i * N + j] = 0.0;

            for (int kk = 0; kk < K; kk += tk) {
                const int ke = imin(kk + tk, K);
                for (int i = ii; i < ie; i++) {
                    double *Ci = C + (size_t)i * N;
                    for (int k = kk; k < ke; k++) {
                        const double a = A[(size_t)i * K + k];
                        const double *Bk = B + (size_t)k * N;
                        for (int j = jj; j < je; j++) Ci[j] += a * Bk[j];
                    }
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr, "usage: %s M K N pr pc b [reps]\n", argv[0]);
        return 1;
    }
    const int M = atoi(argv[1]), K = atoi(argv[2]), N = atoi(argv[3]);
    const int pr = atoi(argv[4]), pc = atoi(argv[5]);
    const int bs = atoi(argv[6]);
    const int reps = (argc > 7) ? atoi(argv[7]) : 10;

    const bool sequential = (pr == 0 && pc == 0);
    if (!sequential && (pr < 1 || pc < 1)) {
        fprintf(stderr, "pr and pc must both be 0 (sequential) or both >= 1\n");
        return 1;
    }
    if (M < 1 || K < 1 || N < 1 || bs < 0 || reps < 1) {
        fprintf(stderr, "invalid arguments\n");
        return 1;
    }
    const int p = sequential ? 1 : pr * pc;

    // Allocation and input generation: outside the timed region.
    std::vector<double> A((size_t)M * K), B((size_t)K * N);
    std::vector<double> C((size_t)M * N), Cref((size_t)M * N);
    gen(A.data(), M, K, 1);
    gen(B.data(), K, N, 2);

    // Reference: sequential UNTILED product (b = 0).
    kernel(A.data(), B.data(), Cref.data(), K, N, 0, M, 0, N, 0);

    const double *Ap = A.data();
    const double *Bp = B.data();
    double *Cp = C.data();

    // Records which thread ran each block, and the team size actually used.
    std::vector<int> owner(p, -1);
    int *owner_p = owner.data();
    int seen = 0;
    int *seen_p = &seen;

    const int rows = sequential ? M : M / pr;
    const int cols = sequential ? N : N / pc;

    std::vector<double> times;
    bool all_ok = true;
    bool mapped_ok = true;

    for (int run = 0; run <= reps; run++) {  // run 0 = untimed warm-up
        std::fill(C.begin(), C.end(), -1.0);  // poison: catches unwritten elements
        std::fill(owner.begin(), owner.end(), -1);
        seen = 0;

        double t0, t1;
        if (sequential) {
            t0 = omp_get_wtime();
            kernel(Ap, Bp, Cp, K, N, 0, M, 0, N, bs);
            t1 = omp_get_wtime();
        } else {
            omp_set_num_threads(p);  // before the timed region
            t0 = omp_get_wtime();
#pragma omp parallel for schedule(static) shared(Ap, Bp, Cp, owner_p, seen_p)
            for (int t = 0; t < p; t++) {
                const int bi = t / pc, bj = t % pc;
                const int r0 = bi * rows;
                const int r1 = (bi == pr - 1) ? M : r0 + rows;
                const int c0 = bj * cols;
                const int c1 = (bj == pc - 1) ? N : c0 + cols;
                owner_p[t] = omp_get_thread_num();
                if (t == 0) *seen_p = omp_get_num_threads();
                kernel(Ap, Bp, Cp, K, N, r0, r1, c0, c1, bs);
            }
            t1 = omp_get_wtime();
        }

        // Exact comparison (entries are small integers, product exact in double).
        if (std::memcmp(C.data(), Cref.data(), C.size() * sizeof(double)) != 0)
            all_ok = false;

        // Iteration t must have been run by thread t, by a team of exactly p.
        if (!sequential) {
            if (seen != p) mapped_ok = false;
            for (int t = 0; t < p; t++)
                if (owner[t] != t) mapped_ok = false;
        }

        if (run > 0) times.push_back(t1 - t0);
    }

    std::vector<double> s = times;
    std::sort(s.begin(), s.end());
    double sum = 0;
    for (double t : times) sum += t;
    const double avg = sum / times.size();
    const size_t n = s.size();
    const double med = (n % 2) ? s[n / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]);
    const double mn = s.front(), mx = s.back();
    const double gflops = 2.0 * M * (double)K * N / med / 1e9;

    printf("# times:");
    for (double t : times) printf(" %.6f", t);
    printf("\n");
    printf("RESULT M=%d K=%d N=%d pr=%d pc=%d p=%d b=%d reps=%d "
           "avg=%.6f median=%.6f min=%.6f max=%.6f gflops=%.3f verified=%s mapped=%s\n",
           M, K, N, pr, pc, p, bs, reps, avg, med, mn, mx, gflops,
           all_ok ? "yes" : "NO", mapped_ok ? "yes" : "NO");

    return (all_ok && mapped_ok) ? 0 : 2;
}