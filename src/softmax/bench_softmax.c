// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <bench_softmax_utils.h>


/** Display the content of a binary32 n-element array */
void array_dump_fp32(float *array, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        printf(" %.3f ", array[i]);
    }
    printf("\n");
}

/** Display the content of a binary64 n-element array */
void array_dump_fp64(double *array, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        printf(" %.3f ", array[i]);
    }
    printf("\n");
}

/** Declaring various softmax implementation benchmarks **/
softmax_bench_result_t softmax_baseline_fp32_bench(float* dst, float* src, double* golden, size_t n);


// Defining a default size fot the inputs and output array
// (can be overloaded during compilation with -DARRAY_SIZE=<value>)
#ifndef ARRAY_SIZE
#define ARRAY_SIZE 16
#endif

float src[ARRAY_SIZE];
float dst[ARRAY_SIZE] = {0.f};
double golden[ARRAY_SIZE] = {0.};

typedef softmax_bench_result_t (softmax_bench_func_t)(float* dst, float* src, double* golden, size_t n);

/** Descriptor structure for softmax benchmark */
typedef struct {
    softmax_bench_func_t* bench;
    char label[100];
} softmax_bench_t;


extern void softmax_baseline_fp32_fp64(double* dst, float* src, size_t n);

int main(void) {
    int i;
    unsigned long start, stop;
    // random initialization of the input arrays
    for (i = 0; i < ARRAY_SIZE; ++i) {
        src[i] = 2.f * rand() / (float) RAND_MAX - 1.f;
    }

    printf("source matrix:\n");
    array_dump_fp32(src, ARRAY_SIZE);

    // computing golden value
    softmax_baseline_fp32_fp64(golden, src, ARRAY_SIZE);
    printf("golden result:\n");
    array_dump_fp64(golden, ARRAY_SIZE);

    softmax_bench_t benchmarks[] = {
        (softmax_bench_t){.bench = softmax_baseline_fp32_bench, .label="baseline n-element softmax"},
    };

    // softmax benchmarks
    for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(softmax_bench_t); benchId++)
    {
        memset(dst, 0, sizeof(dst)); // resetting array in-between experiments
        softmax_bench_result_t bench_result = benchmarks[benchId].bench(dst, src, golden, ARRAY_SIZE);

        printf("--------------------------------------------------------------------------------\n");
        printf("%s result:\n", benchmarks[benchId].label);
        array_dump_fp32(dst, ARRAY_SIZE);

        printf("%s used %d " PERF_METRIC "(s) to evaluate softmax on a %d-element array.\n",
            benchmarks[benchId].label, bench_result.perf_count, ARRAY_SIZE);
        printf("  max absolute error: %.4a\n", bench_result.max_abs_error);
        printf("  max relative error: %.4a\n", bench_result.max_rel_error);
        printf("  error norm 2:       %.4a\n", bench_result.error_norm2);
    }

#if 0
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
#endif

    return 0;
}
