// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <bench_matrix_utils.h>


/** Display the content of a n x matrix on stdout */
void matrix_dump(float *mat, unsigned n)
{
    size_t i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) printf(" %.3f ", mat[i * n + j]);
        printf("\n");
    }
}

/** Declaring various matrix transpose implementations **/


unsigned long matrix_transpose_nxn_bench(float *dst, float *src, size_t n);

unsigned long matrix_transpose_4x4_bench(float *dst, float *src, size_t n); 

unsigned long matrix_transpose_intrinsics_4x4_bench(float* dst, float* src);

unsigned long matrix_transpose_intrinsics_nxn_bench(float *dst, float *src, size_t n);

unsigned long matrix_transpose_intrinsics_loads_nxn_bench(float *dst, float *src, size_t n); 

unsigned long matrix_4x4_transpose_segmented_load_intrinsics_bench(float* dst, float* src);

unsigned long matrix_4x4_transpose_segmented_store_intrinsics_bench(float* dst, float* src);

unsigned long matrix_4x4_transpose_vrgather_bench (float* dst, float* src);

unsigned long matrix_4x4_transpose_vslide_bench(float* dst, float* src); 


// Defining a default size fot the inputs and output array
// (can be overloaded during compilation with -DARRAY_SIZE=<value>)
#ifndef MATRIX_SIZE
#define MATRIX_SIZE 4
#endif

float src[MATRIX_SIZE * MATRIX_SIZE];
float dst[MATRIX_SIZE * MATRIX_SIZE] = {0.f};

typedef unsigned long (matrix_transpose_4x4_bench_func_t)(float* dst, float* src);
typedef unsigned long (matrix_transpose_nxn_bench_func_t)(float* dst, float* src, size_t n);

/** Descriptor structure for 4x4 matrix transpose benchmark */
typedef struct {
    matrix_transpose_4x4_bench_func_t* bench;
    char label[100];
} matrix_4x4_bench_t;

/** Descriptor structure for nxn matrix transpose benchmark */
typedef struct {
    matrix_transpose_nxn_bench_func_t* bench;
    char label[100];
} matrix_nxn_bench_t;


int main(void) {
    int i;
    unsigned long start, stop;
    // random initialization of the input arrays
    for (i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        src[i] = rand() / (float) RAND_MAX;
    }

    printf("source matrix:\n");
    matrix_dump(src, MATRIX_SIZE);

    matrix_4x4_bench_t benchmarks_4x4[] = {
        (matrix_4x4_bench_t){.bench = matrix_transpose_4x4_bench, .label="baseline matrix_transpose 4x4"},
        (matrix_4x4_bench_t){.bench = matrix_transpose_intrinsics_4x4_bench, .label="matrix_transpose intrinsics 4x4"},
        (matrix_4x4_bench_t){.bench = matrix_4x4_transpose_segmented_load_intrinsics_bench, .label="matrix_transpose_segmented_load 4x4"},
        (matrix_4x4_bench_t){.bench = matrix_4x4_transpose_segmented_store_intrinsics_bench, .label="matrix_transpose_segmented_store 4x4"},
        (matrix_4x4_bench_t){.bench = matrix_4x4_transpose_vrgather_bench, .label="matrix_transpose_vrgather 4x4"},
        (matrix_4x4_bench_t){.bench = matrix_4x4_transpose_vslide_bench, .label="matrix_transpose_vslide 4x4"},
    };

    matrix_nxn_bench_t benchmarks_nxn[] = {
        (matrix_nxn_bench_t){.bench = matrix_transpose_nxn_bench, .label="baseline matrix_transpose nxn"},
        (matrix_nxn_bench_t){.bench = matrix_transpose_intrinsics_nxn_bench, .label="matrix_transpose intrinsics nxn"},
        (matrix_nxn_bench_t){.bench = matrix_transpose_intrinsics_loads_nxn_bench, .label="matrix_transpose_loads intrinsics nxn"},
    };

    // 4x4 benchmarks
    for (unsigned benchId=0; benchId < sizeof(benchmarks_4x4) / sizeof(matrix_4x4_bench_t); benchId++)
    {
        memset(dst, 0, sizeof(dst)); // resetting array in-between experiments
        unsigned long perf_count = benchmarks_4x4[benchId].bench(dst, src);

        printf("--------------------------------------------------------------------------------\n");
        printf("%s result:\n", benchmarks_4x4[benchId].label);
        matrix_dump(dst, MATRIX_SIZE);

        printf("%s used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
            benchmarks_4x4[benchId].label, perf_count, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);
    }

    // nxn benchmarks
    for (unsigned benchId=0; benchId < sizeof(benchmarks_nxn) / sizeof(matrix_nxn_bench_t); benchId++)
    {
        memset(dst, 0, sizeof(dst)); // resetting array in-between experiments
        unsigned long perf_count = benchmarks_nxn[benchId].bench(dst, src, MATRIX_SIZE);

        printf("--------------------------------------------------------------------------------\n");
        printf("%s result:\n", benchmarks_nxn[benchId].label);
        matrix_dump(dst, MATRIX_SIZE);

        printf("%s used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
            benchmarks_nxn[benchId].label, perf_count, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);
    }

    return 0;
}
