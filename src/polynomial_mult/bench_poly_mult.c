// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <bench_poly_mult_utils.h>


/** Number of executed test(s) */
#ifndef NUM_TESTS
#define NUM_TESTS 10
#endif


/** Declaring various softmax implementation benchmarks **/
poly_mult_bench_result_t poly_mult_mod_baseline_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_mod_scalar_opt_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_fast_ntt_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_recursive_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_strided_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_strided_barrett_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden); 

poly_mult_bench_result_t poly_mult_ntt_rvv_indexed_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_compressed_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_compressed_barrett_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_indexed_barrett_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

poly_mult_bench_result_t poly_mult_ntt_rvv_fastest_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden); 

typedef poly_mult_bench_result_t (poly_mult_bench_func_t)(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden);

/** Descriptor structure for softmax benchmark */
typedef struct {
    poly_mult_bench_func_t* bench;
    poly_mult_bench_result_t result;
    char label[256];
} poly_mult_bench_t;


// for use as golden implementations
extern void poly_mult_mod_baseline(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo);

/** Display the content of a polynomial on stdout */
void poly_dump(polynomial_t poly)
{
    int i = 0;
    printf(" %d . X^%d ", poly.coeffs[i], i);
    for (i = 1; i <= poly.degree; ++i) {
        printf("+ %d . X^%d ", poly.coeffs[i], i);
    }
    printf("\n");
}

int main(void) {
    int i;
    poly_mult_bench_t benchmarks[] = {
        (poly_mult_bench_t){.bench = poly_mult_mod_baseline_bench,               .label="baseline polynomial multiplication"},
        (poly_mult_bench_t){.bench = poly_mult_ntt_bench,                        .label="slow NTT multiplication"},
        (poly_mult_bench_t){.bench = poly_mult_fast_ntt_bench,                   .label="fast NTT multiplication"},
        (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_bench,                    .label="RVV NTT multiplication"},
        (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_recursive_bench,          .label="RVV NTT multiplication recursive"},
        // (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_indexed_bench,            .label="RVV NTT multiplication split-loops no-recursion indexed-load"},
        // (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_strided_bench,            .label="RVV NTT multiplication split-loops no-recursion strided-load"},
        (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_strided_barrett_bench,    .label="RVV NTT multiplication split-loops no-recursion strided-load with Barrett reduction"},
        // (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_compressed_bench,         .label="RVV NTT multiplication split-loops no-recursion vcompress-based"},
        (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_compressed_barrett_bench, .label="RVV NTT multiplication split-loops no-recursion vcompress-based with Barrett reduction"}, 
        (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_indexed_barrett_bench,    .label="RVV NTT multiplication split-loops no-recursion indexed-load-based with Barrett reduction"}, 
#if LMUL > 1
        (poly_mult_bench_t){.bench = poly_mult_ntt_rvv_fastest_bench, .label="RVV-based ntt-based multiplication (fastest variant [hopefully]"}, 
#endif
    };
    int moduloCoeffs[129] = {0};
    moduloCoeffs[0] = -1;
    moduloCoeffs[128] = 1;
    polynomial_t degree128modulo = {.modulo = 3329, .degree = 128, .coeffSize = 129, .coeffs = moduloCoeffs};

    size_t testSizes[] = {127};//, 16, 255}; // , 16, 17, 32, 33, 128, 129, 511, 512, 1024, 2048};
    polynomial_t modulos[] = {degree128modulo};
    for (size_t testId = 0; testId < sizeof(testSizes) / sizeof(size_t); testId++)
    {
        size_t n = testSizes[testId];
#       ifdef VERBOSE 
        printf("--------------------------------------------------------------------------------\n");
        printf("--------------------------------------------------------------------------------\n");
        printf("Benchmarking polynomial multiplication (with modulo) on %d-degree polynomials.\n", n);
#       endif
        int modulo = 3329; // "small"-Kyber value (8380417 for Dilithium)
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
            poly_mult_mod_baseline(&golden, lhs, rhs, modulos[testId]);

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
#               ifdef VERY_VERBOSE
                printf("running method: %s\n", benchmarks[benchId].label);
#               endif // VERY_VERBOSE

                polynomial_t lhs_copy = copy_poly(lhs);
                polynomial_t rhs_copy = copy_poly(rhs);

                poly_mult_bench_result_t local_result = benchmarks[benchId].bench(&dst, &lhs_copy, &rhs_copy, &modulos[testId], &golden);

                benchmarks[benchId].result = accumulate_bench_result(benchmarks[benchId].result, local_result);

                free(lhs_copy.coeffs);
                free(rhs_copy.coeffs);

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
            printf("%s used %d " PERF_METRIC "(s) to evaluate multiplication on a degree %d polynomial.\n",
                benchmarks[benchId].label, bench_result.perf_count, n);
            printf(" " PERF_METRIC " per multiplication:  %d\n", bench_result.perf_count);
            printf(" " PERF_METRIC " per degree:          %.3f\n", (double) bench_result.perf_count / n);
            printf("  element(s) per " PERF_METRIC ":     %.2e\n", (double) n / bench_result.perf_count);
            printf("  error(s):  %d\n", bench_result.errors);
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
