/* Simple header for ISOC99 output from tce.py
 * Copyright (C) 2013 ET International, Inc.
 * written by Mark Glines
 */
#ifndef _TCE_H
#define _TCE_H

#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <sys/types.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>


/* constants */
#define irrep_e 0
#define irrep_e2 0
#define irrep_f 0
#define irrep_v 0
#define irrep_t 0
#define irrep_t1 0
#define irrep_t2 0
#define irrep_t3 0

/* global GA allocation index things */
extern int *k_spin;
extern int *k_sym;
extern int *k_range;
extern int noab;
extern int nvab;
extern bool restricted, intorb;
extern double FACTORIAL[20];

/* If restricted == true, impose bounds checking on 2d array indices */
extern void tce_restricted_2(/*  inputs */ int  a1b, int  a2b,
                             /* outputs */ int *b1b, int *b2b);
/* If restricted == true, impose bounds checking on 4d array indices */
extern void tce_restricted_4(/*  inputs */ int  a1b, int  a2b, int  a3b, int  a4b,
                             /* outputs */ int *b1b, int *b2b, int *b3b, int *b4b);
extern void tce_sort_2(double *unsorted, double *sorted, int a, int b, int i, int j, double factor);
extern void tce_sort_4(double *unsorted, double *sorted, int a, int b, int c, int d, int i, int j, int k, int l, double factor);
int tce_hash(int *hash, int key);
extern void tce_get_hash_block(double *in, double *out, int dim, int *hash, int key);
extern void tce_add_hash_block(double *out, double *in, int dim, int *hash, int key);
void tce_get_block_ind_i(double *in, double *out, int dim, int w2b, int w1b, int w4b, int w3b);

/* Do a matrix multiply. */
#ifdef USE_MKL_BLAS
#include <mkl_cblas.h>
/* void cblas_dgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                    const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                    const MKL_INT K, const double alpha, const double *A,
                    const MKL_INT lda, const double *B, const MKL_INT ldb,
                    const double beta, double *C, const MKL_INT ldc); */

#define TCE_DGEMM(tA, tB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) cblas_dgemm(CblasColMajor, tA, tB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC)
#elif USE_ATLAS_BLAS
#include <atlas/cblas.h>
/* void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                    const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                    const int K, const double alpha, const double *A,
                    const int lda, const double *B, const int ldb,
                    const double beta, double *C, const int ldc); */
#define TCE_DGEMM(tA, tB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) cblas_dgemm(CblasColMajor, tA, tB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC)
#elif USE_OPENBLAS
#include <cblas.h>
/* void cblas_dgemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                    OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                    OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A,
                    OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb,
                    OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc); */
#define TCE_DGEMM(tA, tB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) cblas_dgemm(CblasColMajor, tA, tB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC)
#else
#warning no BLAS library specified
extern void TCE_DGEMM(char transposeA, char transposeB, int m, int n, int k, double alpha,
        double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);
#endif

#define debug(a, b...) do { printf("%s:%d: " a "\n", __FILE__, __LINE__, ##b); fflush(stdout); } while(0)
#define die(a, b...) do { fprintf(stderr,"FATAL: %s:%d: " a "\n", __FILE__, __LINE__, ##b); fflush(stderr); abort(); } while(0)
#define die_if(a) if((a)) die("failed assert: " #a)

static inline void *_tce_malloc(size_t size) {
    void *rv = calloc(size, 1);
    if(!rv)
        die("calloc of size %zd failed with errno %d (%s)", size, errno, strerror(errno));
    return rv;
}

static inline int *tce_int_malloc(size_t count) {
    int *rv = (int*)_tce_malloc(count * sizeof(int));
    return rv;
}

static inline double *tce_double_malloc(size_t count) {
    double *rv = (double*)_tce_malloc(count * sizeof(double));
    return rv;
}

static inline void tce_free(void *ptr) {
    free(ptr);
}

void write_tensor(char *filename, double *data, size_t len);
void dump_tensor(double *data, int *hash, int dims, int *sizes, int *starts);
void dump_block(double *data, int len);

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

#endif /* _TCE_H */
