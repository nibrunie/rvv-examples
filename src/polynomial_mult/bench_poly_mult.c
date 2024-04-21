// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <bench_poly_mult_utils.h>


/** Declaring various softmax implementation benchmarks **/
poly_mult_bench_result_t poly_mult_baseline_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* golden);

typedef poly_mult_bench_result_t (poly_mult_bench_func_t)(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* golden);

/** Descriptor structure for softmax benchmark */
typedef struct {
    poly_mult_bench_func_t* bench;
    poly_mult_bench_result_t result;
    char label[100];
} poly_mult_bench_t;


// for use as golden implementation
extern void poly_mult_baseline(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs);

/** Display the content of a polynomial on stdout */
void poly_dump(polynomial_t poly)
{
    int i;
    for (i = 0; i <= poly.degree; ++i) {
        printf(" %d . X^%d ", poly.coeffs[i], i);
    }
    printf("\n");
}

#ifndef NUM_TESTS
#define NUM_TESTS 100
#endif


int main(void) {
    int i;
    poly_mult_bench_t benchmarks[] = {
        (poly_mult_bench_t){.bench = poly_mult_baseline_bench, .label="baseline polynomial multiplication"},
    };

    size_t testSizes[] = {4}; // , 16, 17, 32, 33, 128, 129, 511, 512, 1024, 2048};
    for (size_t testId = 0; testId < sizeof(testSizes) / sizeof(size_t); testId++)
    {
        size_t n = testSizes[testId];
#       ifdef VERBOSE 
        printf("--------------------------------------------------------------------------------\n");
        printf("--------------------------------------------------------------------------------\n");
        printf("Benchmarking polynomial on %d-degree polynomials.\n", n);
#       endif
        int modulo = 8319;
        polynomial_t lhs = allocate_poly(n, modulo);
        polynomial_t rhs = allocate_poly(n, modulo);
        polynomial_t golden = allocate_poly(2*n, modulo);
        polynomial_t dst = allocate_poly(2*n, modulo);
        assert(lhs.coeffs);
        assert(rhs.coeffs);
        assert(golden.coeffs);
        assert(dst.coeffs);

        // reset benchmark results
        for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(poly_mult_bench_t); benchId++)
        {
            benchmarks[benchId].result.errors = 0;
            benchmarks[benchId].result.perf_count = 0;
        }

        int j;

        for (j = 0; j < NUM_TESTS; ++j) {
            randomize_poly(rhs);
            randomize_poly(lhs);

            // computing golden value
            poly_mult_baseline(&golden, lhs, rhs);

#           ifdef VERY_VERBOSE
            printf("lhs polynomial :\n");
            poly_dump(lhs);
            printf("rhs polynomial :\n");
            poly_dump(rhs);
            printf("golden result:\n");
            poly_dump(golden);
#           endif // VERY_VERBOSE

            // softmax benchmarks. iterating over all existing implementation for this given input set
            for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(poly_mult_bench_t); benchId++)
            {
                memset(dst.coeffs, 0, sizeof(uint32_t) * dst.degree); // resetting array in-between experiments

                poly_mult_bench_result_t local_result = benchmarks[benchId].bench(&dst, &lhs, &rhs, &golden);

                benchmarks[benchId].result = accumulate_bench_result(benchmarks[benchId].result, local_result);

#               ifdef VERY_VERBOSE
                printf("%s result:\n", benchmarks[benchId].label);
                poly_dump(dst);
#               endif // VERY_VERBOSE

            }
        }

        // display results
        for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(poly_mult_bench_t); benchId++)
        {
            poly_mult_bench_result_t bench_result = benchmarks[benchId].result;
            bench_result.perf_count = bench_result.perf_count / NUM_TESTS;


#           ifdef VERBOSE 
            printf("--------------------------------------------------------------------------------\n");
            printf("%s used %d " PERF_METRIC "(s) to evaluate softmax on a %d-element array.\n",
                benchmarks[benchId].label, bench_result.perf_count, n);
            printf(" " PERF_METRIC " per elements:    %.3f\n", (double) bench_result.perf_count / n);
            printf("  element(s) per " PERF_METRIC ": %.3f\n", (double) n / bench_result.perf_count);
            printf("  errors:  %.4a\n", bench_result.errors);
#           else
            // condensed display
            printf("%s, %d, %d, %d\n", 
                   benchmarks[benchId].label, n, bench_result.perf_count,
                   bench_result.errors);
#           endif
        }

        free(lhs.coeffs);
        free(rhs.coeffs);
        free(dst.coeffs);
        free(golden.coeffs);
    }


    return 0;
}
