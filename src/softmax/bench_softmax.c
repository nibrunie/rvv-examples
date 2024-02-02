// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <bench_matrix_utils.h>


/** Display the content of a n x matrix on stdout */
void array_dump(float *array, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        printf(" %.3f ", array[i]);
    }
    printf("\n");
}

/** Declaring various matrix transpose implementations **/
void softmax_baseline_fp32(float*, float*, size_t);


// Defining a default size fot the inputs and output array
// (can be overloaded during compilation with -DARRAY_SIZE=<value>)
#ifndef ARRAY_SIZE
#define ARRAY_SIZE 4
#endif

float src[ARRAY_SIZE];
float dst[ARRAY_SIZE] = {0.f};

typedef struct {
    unsigned long cycles;
    unsigned long instret;
    double accuracy;
} softmax_bench_result_t;

typedef softmax_bench_result_t (softmax_bench_func_t)(float* dst, float* src, size_t n);

/** Descriptor structure for softmax benchmark */
typedef struct {
    softmax_bench_func_t* bench;
    char label[100];
} softmax_bench_t;



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

        printf("%s used %d " PERF_METRIC "(s) to tranpose %dx%d=%d element(s).\n",
            benchmarks_4x4[benchId].label, perf_count, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);
    }

    // nxn benchmarks on 4x4 tranpose
    for (unsigned benchId=0; benchId < sizeof(benchmarks_nxn) / sizeof(matrix_nxn_bench_t); benchId++)
    {
        memset(dst, 0, sizeof(dst)); // resetting array in-between experiments
        unsigned long perf_count = benchmarks_nxn[benchId].bench(dst, src, MATRIX_SIZE);

        printf("--------------------------------------------------------------------------------\n");
        printf("%s result:\n", benchmarks_nxn[benchId].label);
        matrix_dump(dst, MATRIX_SIZE);

        printf("%s used %d " PERF_METRIC "(s) to tranpose %dx%d=%d element(s).\n",
            benchmarks_nxn[benchId].label, perf_count, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);
    }
    size_t testSizes[] = {4, 16, 128, 512, 17, 129, 511};
    for (size_t testId = 0; testId < sizeof(testSizes) / sizeof(size_t); testId++)
    {
        size_t n = testSizes[testId];
        float* matrixIn = malloc(n * n * sizeof(float));
        assert(matrixIn);

        float* matrixOut = malloc(n * n * sizeof(float));
        assert(matrixOut);

        float* matrixRef = malloc(n * n * sizeof(float));
        assert(matrixRef);


        for (i = 0; i < n * n; ++i) {
            matrixIn[i] = rand() / (float) RAND_MAX;
        }

        // computing reference
        matrix_transpose_nxn_bench(matrixRef, matrixIn, n);

        for (unsigned benchId=0; benchId < sizeof(benchmarks_nxn) / sizeof(matrix_nxn_bench_t); benchId++)
        {
            memset(matrixOut, 0, n * n * sizeof(float)); // resetting array in-between experiments
            unsigned long perf_count = benchmarks_nxn[benchId].bench(matrixOut, matrixIn, n);

            // comparing with golden (comparison mute for matrix_transpose_nxn_bench which is used as golden
            assert(!memcmp(matrixRef, matrixOut, n * n * sizeof(float)));

            printf("--------------------------------------------------------------------------------\n");

            printf("%s used %d " PERF_METRIC "(s) to tranpose %dx%d=%d element(s).\n",
                benchmarks_nxn[benchId].label, perf_count, n, n, n * n);
        }

        free(matrixIn);
        free(matrixOut);
        free(matrixRef);

    }

    return 0;
}
