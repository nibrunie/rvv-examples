// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <bench_softmax_utils.h>


/** Declaring various softmax implementation benchmarks **/
softmax_bench_result_t softmax_baseline_fp32_bench(float* dst, float* src, double* golden, size_t n);

softmax_bench_result_t softmax_rvv_norm_fp32_bench(float* dst, float* src, double* golden, size_t n); 

softmax_bench_result_t softmax_rvv_fp32_bench(float* dst, float* src, double* golden, size_t n); 

softmax_bench_result_t softmax_accurate_rvv_fp32_bench(float* dst, float* src, double* golden, size_t n);

typedef softmax_bench_result_t (softmax_bench_func_t)(float* dst, float* src, double* golden, size_t n);

/** Descriptor structure for softmax benchmark */
typedef struct {
    softmax_bench_func_t* bench;
    char label[100];
} softmax_bench_t;


extern void softmax_baseline_fp32_fp64(double* dst, float* src, size_t n);


int main(void) {
    int i;
    softmax_bench_t benchmarks[] = {
        (softmax_bench_t){.bench = softmax_baseline_fp32_bench,     .label="baseline n-element softmax"},
        (softmax_bench_t){.bench = softmax_rvv_norm_fp32_bench,     .label="rvv-based (norm only) n-element softmax"},
        (softmax_bench_t){.bench = softmax_rvv_fp32_bench,          .label="rvv-based n-element softmax"},
        (softmax_bench_t){.bench = softmax_accurate_rvv_fp32_bench, .label="rvv-based n-element accurate softmax"},
    };

    size_t testSizes[] = {4, 16, 17, 128, 129, 511, 512, 1024, 2048};
    for (size_t testId = 0; testId < sizeof(testSizes) / sizeof(size_t); testId++)
    {
        size_t n = testSizes[testId];
#       ifdef VERBOSE 
        printf("--------------------------------------------------------------------------------\n");
        printf("--------------------------------------------------------------------------------\n");
        printf("Benchmarking softmax on a %d-element array.\n", n);
#       endif
        float* src = malloc(n * sizeof(float));
        float* dst = malloc(n * sizeof(float));
        double* golden = malloc(n * sizeof(double));
        assert(src);
        assert(dst);
        assert(golden);
        // random initialization of the input arrays
        for (i = 0; i < n; ++i) {
            src[i] = 2.f * rand() / (float) RAND_MAX - 1.f;
        }

        // printf("source matrix:\n");
        // array_dump_fp32(src, n);

        // computing golden value
        softmax_baseline_fp32_fp64(golden, src, n);
        // printf("golden result:\n");
        // array_dump_fp64(golden, n);


        // softmax benchmarks
        for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(softmax_bench_t); benchId++)
        {
            memset(dst, 0, sizeof(dst)); // resetting array in-between experiments
            softmax_bench_result_t bench_result = benchmarks[benchId].bench(dst, src, golden, n);


#           ifdef VERBOSE 
            printf("--------------------------------------------------------------------------------\n");
            // printf("%s result:\n", benchmarks[benchId].label);
            // array_dump_fp32(dst, n);
            printf("%s used %d " PERF_METRIC "(s) to evaluate softmax on a %d-element array.\n",
                benchmarks[benchId].label, bench_result.perf_count, n);
            printf(" " PERF_METRIC " per elements:    %.3f\n", (double) bench_result.perf_count / n);
            printf("  element(s) per " PERF_METRIC ": %.3f\n", (double) n / bench_result.perf_count);
            printf("  max absolute error: %.4a\n", bench_result.max_abs_error);
            printf("  max relative error: %.4a\n", bench_result.max_rel_error);
            printf("  error norm 2:       %.4a\n", bench_result.error_norm2);
#           else
            // condensed display
            printf("%s, %d, %d, %.3e, %.3e, %.3e\n", 
                   benchmarks[benchId].label, n, bench_result.perf_count,
                   bench_result.max_abs_error, bench_result.max_rel_error, bench_result.error_norm2);
#           endif
        }

        free(src);
        free(dst);
        free(golden);
    }


    return 0;
}
