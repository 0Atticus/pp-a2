/* assignment2.h: starter kit for ITCS 4145 Assignment 2 (Parts A and B).
 *
 * Include this header before any other header, since it sets the feature
 * macro that makes clock_gettime visible. Build pthreads programs with
 * gcc -O2 -pthread and OpenMP programs with gcc -O2 -fopenmp. The OpenMP
 * programs time with omp_get_wtime, not with now_seconds.
 *
 * Contents:
 *   now_seconds()       monotonic wall-clock time in seconds (Part A)
 *   generate_matrix()   row-major dense matrix of small random integers
 *
 * All functions are static inline, so each .c file that includes this
 * header gets its own copy, no separate object file is needed, and unused
 * functions produce no warning under -Wall.
 */
#ifndef ASSIGNMENT2_H
#define ASSIGNMENT2_H

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* ---------- Timer -------------------------------------------------------- */

/* Seconds since an arbitrary fixed point; differences are elapsed wall time. */
static inline double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ---------- Random numbers ---------------------------------------------- */

/* xorshift64*: the same seed gives the same sequence on every machine. */
typedef struct { uint64_t state; } rng_t;

static inline void rng_seed(rng_t *r, uint64_t seed) {
    r->state = seed * 0x9E3779B97F4A7C15ULL + 1;
} 

static inline uint64_t rng_next(rng_t *r) {
    r->state ^= r->state >> 12;
    r->state ^= r->state << 25;
    r->state ^= r->state >> 27;
    return r->state * 2685821657736338717ULL;
}

/* Uniform double in [0, 1). */
static inline double rng_unit(rng_t *r) {
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}

/* ---------- Dense matrices ------------------------------------------------ */

/* A rows x cols matrix in row-major order, element (i, j) at a[i*cols + j],
 * filled with random integers from 1 to 9 stored as doubles. Products and
 * sums of such entries are exact in double for the sizes of this
 * assignment, so results can be compared for exact equality.
 * Exits with a message if allocation fails. */
static inline double *generate_matrix(int rows, int cols, uint64_t seed) {
    double *a = (double *)malloc((size_t)rows * (size_t)cols * sizeof(double));
    if (a == NULL) {
        fprintf(stderr, "generate_matrix: out of memory\n");
        exit(1);
    }
    rng_t r;
    rng_seed(&r, seed);
    size_t count = (size_t)rows * (size_t)cols;
    for (size_t i = 0; i < count; i++)
        a[i] = (double)(1 + (int)(rng_next(&r) % 9));
    return a;
}

#endif /* ASSIGNMENT2_H */
