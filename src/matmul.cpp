// ITCS 4145 Assignment 2, Part A: blocked matrix multiplication with pthreads.
//
// Usage: ./matmul M K N pr pc b [reps]
//   pr = pc = 0 : sequential implementation (no threads created)
//   pr x pc     : thread grid (pr x 1 is the 1D row partition)
//   b           : cache tile size, 0 = no tiling
//   reps        : timed repetitions (default 10), after one untimed warm-up
//
// Build: g++ -O2 -pthread matmul.cpp -o matmul

#include <pthread.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "assignment2.h"

// ---------------------------------------------------------------------------
// ADJUST THESE TWO WRAPPERS to the real names/signatures in assignment2.h.
// ---------------------------------------------------------------------------
static void gen(int rows, int cols, unsigned seed) {
    generate_matrix(rows, cols, seed);  // <-- check signature in header
}
static double wall_time() {
    return now_seconds();                  // <-- check name in header
}
// ---------------------------------------------------------------------------

static inline int imin(int a, int b) { return a < b ? a : b; }

// Tiled kernel on block C[r0,r1) x [c0,c1). b == 0 means one tile (no tiling).
// Loop order inside a tile: i, k, j.
static void kernel(const double *A, const double *B, double *C,
                   int M, int K, int N,
                   int r0, int r1, int c0, int c1, int b) {
    (void)M;
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

struct ThreadArgs {
    const double *A;
    const double *B;
    double *C;
    int M, K, N;
    int r0, r1, c0, c1;
    int b;
};

static void *worker(void *p) {
    const ThreadArgs *t = static_cast<const ThreadArgs *>(p);
    kernel(t->A, t->B, t->C, t->M, t->K, t->N, t->r0, t->r1, t->c0, t->c1, t->b);
    return nullptr;
}

// Fill argument structs (done OUTSIDE the timed region).
static void make_args(std::vector<ThreadArgs> &args, const double *A, const double *B,
                      double *C, int M, int K, int N, int pr, int pc, int b) {
    const int rows = M / pr, cols = N / pc;
    for (int rank = 0; rank < pr * pc; rank++) {
        const int bi = rank / pc, bj = rank % pc;
        ThreadArgs &t = args[rank];
        t.A = A; t.B = B; t.C = C;
        t.M = M; t.K = K; t.N = N;
        t.r0 = bi * rows;
        t.r1 = (bi == pr - 1) ? M : t.r0 + rows;
        t.c0 = bj * cols;
        t.c1 = (bj == pc - 1) ? N : t.c0 + cols;
        t.b = b;
    }
}

int main(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr, "usage: %s M K N pr pc b [reps]\n", argv[0]);
        return 1;
    }
    const int M = atoi(argv[1]), K = atoi(argv[2]), N = atoi(argv[3]);
    const int pr = atoi(argv[4]), pc = atoi(argv[5]);
    const int b = atoi(argv[6]);
    const int reps = (argc > 7) ? atoi(argv[7]) : 10;

    const bool sequential = (pr == 0 && pc == 0);
    if (!sequential && (pr < 1 || pc < 1)) {
        fprintf(stderr, "pr and pc must both be 0 (sequential) or both >= 1\n");
        return 1;
    }
    if (M < 1 || K < 1 || N < 1 || b < 0 || reps < 1) {
        fprintf(stderr, "invalid arguments\n");
        return 1;
    }
    const int p = sequential ? 1 : pr * pc;

    // Allocation and input generation: outside the timed region.
    std::vector<double> A((size_t)M * K), B((size_t)K * N);
    std::vector<double> C((size_t)M * N), Cref((size_t)M * N);
    gen(M, K, 1);
    gen(K, N, 2);

    // Reference: sequential UNTILED product (b = 0).
    kernel(A.data(), B.data(), Cref.data(), M, K, N, 0, M, 0, N, 0);

    std::vector<pthread_t> tid(p);
    std::vector<ThreadArgs> args(sequential ? 0 : p);
    if (!sequential)
        make_args(args, A.data(), B.data(), C.data(), M, K, N, pr, pc, b);

    std::vector<double> times;
    bool all_ok = true;

    for (int run = 0; run <= reps; run++) {  // run 0 = untimed warm-up
        // Poison C so unwritten elements are caught by verification.
        std::fill(C.begin(), C.end(), -1.0);

        double t0 = wall_time();
        if (sequential) {
            kernel(A.data(), B.data(), C.data(), M, K, N, 0, M, 0, N, b);
        } else {
            for (int r = 0; r < p; r++)
                if (pthread_create(&tid[r], nullptr, worker, &args[r]) != 0) {
                    fprintf(stderr, "pthread_create failed (rank %d)\n", r);
                    return 1;
                }
            for (int r = 0; r < p; r++) pthread_join(tid[r], nullptr);
        }
        double t1 = wall_time();

        // Exact comparison (entries are small integers, product exact in double).
        if (std::memcmp(C.data(), Cref.data(), C.size() * sizeof(double)) != 0)
            all_ok = false;

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
           "avg=%.6f median=%.6f min=%.6f max=%.6f gflops=%.3f verified=%s\n",
           M, K, N, pr, pc, p, b, reps, avg, med, mn, mx, gflops,
           all_ok ? "yes" : "NO");

    return all_ok ? 0 : 2;
}