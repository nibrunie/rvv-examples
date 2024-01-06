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


void matrix_transpose(float *dst, float *src, size_t n);

void matrix_transpose_intrinsics_4x4(float *dst, float *src);

void matrix_transpose_intrinsics(float *dst, float *src, size_t n);

void matrix_transpose_intrinsics_loads(float *dst, float *src, size_t n); 

void matrix_4x4_transpose_segmented_load_intrinsics(float* dst, float* src);

void matrix_4x4_transpose_segmented_store_intrinsics(float* dst, float* src);

unsigned long matrix_4x4_transpose_vrgather_bench (float* dst, float* src);

unsigned long matrix_4x4_transpose_vslide_bench(float* dst, float* src); 


// Defining a default size fot the inputs and output array
// (can be overloaded during compilation with -DARRAY_SIZE=<value>)
#ifndef MATRIX_SIZE
#define MATRIX_SIZE 4
#endif

float src[MATRIX_SIZE * MATRIX_SIZE];
float dst[MATRIX_SIZE * MATRIX_SIZE] = {0.f};


int main(void) {
    int i;
    unsigned long start, stop;
    // random initialization of the input arrays
    for (i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        src[i] = rand() / (float) RAND_MAX;
    }

    printf("source matrix:\n");
    matrix_dump(src, MATRIX_SIZE);

    start = read_perf_counter();
    matrix_transpose(dst, src, MATRIX_SIZE);
    stop = read_perf_counter();

    printf("matrix_transpose result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("baseline matrix_transpose used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    memset(dst, 0, sizeof(dst)); // resetting array in-between experiments

    start = read_perf_counter();
    matrix_transpose_intrinsics_4x4(dst, src);
    stop = read_perf_counter();

    printf("matrix_transpose_intrinsics_4x4 result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_transpose_intrinsics_4x4 used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    memset(dst, 0, sizeof(dst)); // resetting array in-between experiments

    start = read_perf_counter();
    matrix_transpose_intrinsics(dst, src, MATRIX_SIZE);
    stop = read_perf_counter();

    printf("matrix_transpose_intrinsics result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_transpose_intrinsics used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    memset(dst, 0, sizeof(dst)); // resetting array in-between experiments

    start = read_perf_counter();
    matrix_transpose_intrinsics_loads(dst, src, MATRIX_SIZE);
    stop = read_perf_counter();

    printf("matrix_transpose_intrinsics_loads result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_transpose_intrinsics_loads used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    memset(dst, 0, sizeof(dst)); // resetting array in-between experiments

    start = read_perf_counter();
    matrix_4x4_transpose_segmented_load_intrinsics(dst, src);
    stop = read_perf_counter();

    printf("matrix_4x4_transpose_segmented_load_intrinscs result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_4x4_transpose_segmented_load_intrinscs used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    memset(dst, 0, sizeof(dst)); // resetting array in-between experiments

    start = read_perf_counter();
    matrix_4x4_transpose_segmented_store_intrinsics(dst, src);
    stop = read_perf_counter();

    printf("matrix_4x4_transpose_segmented_store_intrinscs result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_4x4_transpose_segmented_store_intrinscs used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    memset(dst, 0, sizeof(dst)); // resetting array in-between experiments

    start = 0;
    stop = matrix_4x4_transpose_vrgather_bench(dst, src);

    printf("matrix_4x4_transpose_segmented_store_intrinscs result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_4x4_transpose_vrgather used %d instruction(s) to tranpose %dx%d=%d element(s) (in registers).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    memset(dst, 0, sizeof(dst)); // resetting array in-between experiments

    start = 0;
    stop = matrix_4x4_transpose_vslide_bench(dst, src);

    printf("matrix_4x4_transpose_segmented_store_intrinscs result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_4x4_transpose_vslide used %d instruction(s) to tranpose %dx%d=%d element(s) (in registers).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);
    return 0;
}
