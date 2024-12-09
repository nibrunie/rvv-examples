// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <bench_utils.h>
#include "ubench_utils.h"


/** Number of executed test(s) */
#ifndef NUM_TESTS
#define NUM_TESTS 10
#endif

/** Declaring various softmax implementation benchmarks **/
ubench_result_t baseline_bench(size_t n);


ubench_result_t baseline_bench(size_t n) {
    long start = read_perf_counter();
    for (int i = 0; i < n / 16; i++) {
        asm volatile(
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            "nop\n"
            :
            :
            :
        );
        // do nothing
    }
    long stop = read_perf_counter();
    return (ubench_result_t){
        .perf_count = (stop - start),
        .elt_per_op = 1,
        .errors = 0
    };
};

BENCH_0OP_INSN(nop)
BENCH_LAT_2OP_INSN(add)
BENCH_LAT_2OP_INSN(div)
BENCH_LAT_2OP_INSN(divu)
BENCH_LAT_2OP_INSN(divw)
BENCH_LAT_2OP_INSN(divuw)
BENCH_LAT_2OP_INSN(rem)
BENCH_LAT_2OP_INSN(remu)
BENCH_LAT_2OP_INSN(remw)
BENCH_LAT_2OP_INSN(remuw)
BENCH_LAT_2OP_INSN(mul)
BENCH_LAT_2OP_INSN(mulh)
BENCH_LAT_2OP_INSN(mulhsu)
BENCH_LAT_2OP_INSN(mulhu)
BENCH_LAT_2OP_INSN(mulw)

BENCH_THROUGHPUT_2OP_INSN(add)
BENCH_THROUGHPUT_2OP_INSN(div)
BENCH_THROUGHPUT_2OP_INSN(divu)
BENCH_THROUGHPUT_2OP_INSN(divw)
BENCH_THROUGHPUT_2OP_INSN(divuw)
BENCH_THROUGHPUT_2OP_INSN(rem)
BENCH_THROUGHPUT_2OP_INSN(remu)
BENCH_THROUGHPUT_2OP_INSN(remw)
BENCH_THROUGHPUT_2OP_INSN(remuw)
BENCH_THROUGHPUT_2OP_INSN(mul)
BENCH_THROUGHPUT_2OP_INSN(mulh)
BENCH_THROUGHPUT_2OP_INSN(mulhsu)
BENCH_THROUGHPUT_2OP_INSN(mulhu)
BENCH_THROUGHPUT_2OP_INSN(mulw)


BENCH_LAT_2OP_FPD_INSN(fadd)
BENCH_LAT_2OP_FPD_INSN(fsub)
BENCH_LAT_2OP_FPD_INSN(fdiv)
BENCH_LAT_2OP_FPD_INSN(fmul)

BENCH_THROUGHPUT_2OP_FPD_INSN(fadd)
BENCH_THROUGHPUT_2OP_FPD_INSN(fsub)
BENCH_THROUGHPUT_2OP_FPD_INSN(fdiv)
BENCH_THROUGHPUT_2OP_FPD_INSN(fmul)

BENCH_LAT_2OP_VEC_INSN(vadd, 1, 32)
BENCH_LAT_2OP_VEC_INSN(vadd, 2, 32)
BENCH_LAT_2OP_VEC_INSN(vadd, 4, 32)
BENCH_LAT_2OP_VEC_INSN(vadd, 8, 32)

BENCH_LAT_2OP_VEC_INSN(vfadd, 1, 32)
BENCH_LAT_2OP_VEC_INSN(vfadd, 2, 32)
BENCH_LAT_2OP_VEC_INSN(vfadd, 4, 32)
BENCH_LAT_2OP_VEC_INSN(vfadd, 8, 32)

BENCH_LAT_2OP_VEC_INSN(vfmul, 1, 32)
BENCH_LAT_2OP_VEC_INSN(vfmul, 2, 32)
BENCH_LAT_2OP_VEC_INSN(vfmul, 4, 32)
BENCH_LAT_2OP_VEC_INSN(vfmul, 8, 32)

BENCH_LAT_2OP_VEC_INSN(vfdiv, 1, 32)
BENCH_LAT_2OP_VEC_INSN(vfdiv, 2, 32)
BENCH_LAT_2OP_VEC_INSN(vfdiv, 4, 32)
BENCH_LAT_2OP_VEC_INSN(vfdiv, 8, 32)

BENCH_LAT_2OP_VEC_INSN(vadd, 1, 64)
BENCH_LAT_2OP_VEC_INSN(vadd, 2, 64)
BENCH_LAT_2OP_VEC_INSN(vadd, 4, 64)
BENCH_LAT_2OP_VEC_INSN(vadd, 8, 64)

BENCH_LAT_2OP_VEC_INSN(vfadd, 1, 64)
BENCH_LAT_2OP_VEC_INSN(vfadd, 2, 64)
BENCH_LAT_2OP_VEC_INSN(vfadd, 4, 64)
BENCH_LAT_2OP_VEC_INSN(vfadd, 8, 64)

BENCH_LAT_2OP_VEC_INSN(vfmul, 1, 64)
BENCH_LAT_2OP_VEC_INSN(vfmul, 2, 64)
BENCH_LAT_2OP_VEC_INSN(vfmul, 4, 64)
BENCH_LAT_2OP_VEC_INSN(vfmul, 8, 64)

BENCH_LAT_2OP_VEC_INSN(vfdiv, 1, 64)
BENCH_LAT_2OP_VEC_INSN(vfdiv, 2, 64)
BENCH_LAT_2OP_VEC_INSN(vfdiv, 4, 64)
BENCH_LAT_2OP_VEC_INSN(vfdiv, 8, 64)

int main(void) {
    int i;
    ubench_t benchmarks[] = {
        (ubench_t){.bench = baseline_bench,  .label="baseline benchmark"},
        BENCH_LAT_INSN_TC(nop),
        BENCH_LAT_INSN_TC(add),
        BENCH_LAT_INSN_TC(div),
        BENCH_LAT_INSN_TC(divu),
        BENCH_LAT_INSN_TC(divw),
        BENCH_LAT_INSN_TC(divuw),
        BENCH_LAT_INSN_TC(rem),
        BENCH_LAT_INSN_TC(remu),
        BENCH_LAT_INSN_TC(remw),
        BENCH_LAT_INSN_TC(remuw),
        BENCH_LAT_INSN_TC(mul),
        BENCH_LAT_INSN_TC(mulh),
        BENCH_LAT_INSN_TC(mulhsu),
        BENCH_LAT_INSN_TC(mulhu),
        BENCH_LAT_INSN_TC(mulw),

        BENCH_LAT_INSN_TC(fadd),
        BENCH_LAT_INSN_TC(fsub),
        BENCH_LAT_INSN_TC(fdiv),
        BENCH_LAT_INSN_TC(fmul),

        BENCH_LAT_VEC_INSN_TC(vadd, 1, 32),
        BENCH_LAT_VEC_INSN_TC(vadd, 2, 32),
        BENCH_LAT_VEC_INSN_TC(vadd, 4, 32),
        BENCH_LAT_VEC_INSN_TC(vadd, 8, 32),

        BENCH_LAT_VEC_INSN_TC(vfadd, 1, 32),
        BENCH_LAT_VEC_INSN_TC(vfadd, 2, 32),
        BENCH_LAT_VEC_INSN_TC(vfadd, 4, 32),
        BENCH_LAT_VEC_INSN_TC(vfadd, 8, 32),

        BENCH_LAT_VEC_INSN_TC(vfmul, 1, 32),
        BENCH_LAT_VEC_INSN_TC(vfmul, 2, 32),
        BENCH_LAT_VEC_INSN_TC(vfmul, 4, 32),
        BENCH_LAT_VEC_INSN_TC(vfmul, 8, 32),

        BENCH_LAT_VEC_INSN_TC(vfdiv, 1, 32),
        BENCH_LAT_VEC_INSN_TC(vfdiv, 2, 32),
        BENCH_LAT_VEC_INSN_TC(vfdiv, 4, 32),
        BENCH_LAT_VEC_INSN_TC(vfdiv, 8, 32),

        BENCH_LAT_VEC_INSN_TC(vadd, 1, 64),
        BENCH_LAT_VEC_INSN_TC(vadd, 2, 64),
        BENCH_LAT_VEC_INSN_TC(vadd, 4, 64),
        BENCH_LAT_VEC_INSN_TC(vadd, 8, 64),

        BENCH_LAT_VEC_INSN_TC(vfadd, 1, 64),
        BENCH_LAT_VEC_INSN_TC(vfadd, 2, 64),
        BENCH_LAT_VEC_INSN_TC(vfadd, 4, 64),
        BENCH_LAT_VEC_INSN_TC(vfadd, 8, 64),

        BENCH_LAT_VEC_INSN_TC(vfmul, 1, 64),
        BENCH_LAT_VEC_INSN_TC(vfmul, 2, 64),
        BENCH_LAT_VEC_INSN_TC(vfmul, 4, 64),
        BENCH_LAT_VEC_INSN_TC(vfmul, 8, 64),

        BENCH_LAT_VEC_INSN_TC(vfdiv, 1, 64),
        BENCH_LAT_VEC_INSN_TC(vfdiv, 2, 64),
        BENCH_LAT_VEC_INSN_TC(vfdiv, 4, 64),
        BENCH_LAT_VEC_INSN_TC(vfdiv, 8, 64),

        BENCH_THROUGHPUT_INSN_TC(add),
        BENCH_THROUGHPUT_INSN_TC(div),
        BENCH_THROUGHPUT_INSN_TC(divu),
        BENCH_THROUGHPUT_INSN_TC(divw),
        BENCH_THROUGHPUT_INSN_TC(divuw),
        BENCH_THROUGHPUT_INSN_TC(rem),
        BENCH_THROUGHPUT_INSN_TC(remu),
        BENCH_THROUGHPUT_INSN_TC(remw),
        BENCH_THROUGHPUT_INSN_TC(remuw),
        BENCH_THROUGHPUT_INSN_TC(mul),
        BENCH_THROUGHPUT_INSN_TC(mulh),
        BENCH_THROUGHPUT_INSN_TC(mulhsu),
        BENCH_THROUGHPUT_INSN_TC(mulhu),
        BENCH_THROUGHPUT_INSN_TC(mulw),

        BENCH_THROUGHPUT_INSN_TC(fadd),
        BENCH_THROUGHPUT_INSN_TC(fsub),
        BENCH_THROUGHPUT_INSN_TC(fdiv),
        BENCH_THROUGHPUT_INSN_TC(fmul),
    };
#ifndef VERBOSE
    // condensed display
    printf("label, test-size, counter-value, metric-per-element, element-per-metric, error(s)\n");
#endif

    size_t testSizes[] = {65536};
    for (size_t testId = 0; testId < sizeof(testSizes) / sizeof(size_t); testId++)
    {
        size_t n = testSizes[testId];
#       ifdef VERBOSE 
        printf("--------------------------------------------------------------------------------\n");
        printf("--------------------------------------------------------------------------------\n");
        printf("Benchmarking with test size %d.\n", n);
#       endif

        // reset benchmark results
        for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(ubench_t); benchId++)
        {
            benchmarks[benchId].result.errors = 0;
            benchmarks[benchId].result.perf_count = 0;
        }

        for (int j = 0; j < NUM_TESTS; ++j) {
            // softmax benchmarks. iterating over all existing implementation for this given input set
            for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(ubench_t); benchId++)
            {
#               ifdef VERY_VERBOSE
                printf("running method: %s\n", benchmarks[benchId].label);
#               endif // VERY_VERBOSE

                ubench_result_t local_result = benchmarks[benchId].bench(n);
                benchmarks[benchId].result = accumulate_ubench_result(benchmarks[benchId].result, local_result);
            }
        }

        // display results
        for (unsigned benchId=0; benchId < sizeof(benchmarks) / sizeof(ubench_t); benchId++)
        {
            ubench_result_t bench_result = benchmarks[benchId].result;
            bench_result.perf_count = bench_result.perf_count / NUM_TESTS;


#           ifdef VERBOSE 
            printf("--------------------------------------------------------------------------------\n");
            printf("%s used %d " PERF_METRIC "(s) to evaluate multiplication on a degree %d polynomial.\n",
                benchmarks[benchId].label, bench_result.perf_count, n);
            printf("  " PERF_METRIC " per run:        %d\n", bench_result.perf_count);
            printf("  " PERF_METRIC " per element:    %.3f\n", (double) bench_result.perf_count / (n * bench_result.elt_per_op));
            printf("  element(s) per " PERF_METRIC ": %.2e\n", (double) (n * bench_result.elt_per_op) / bench_result.perf_count);
            printf("  error(s):  %d\n", bench_result.errors);
#           else
            // condensed display
            printf("%s, %d, %d, %.3f, %.2f, %d\n", 
                   benchmarks[benchId].label, n, bench_result.perf_count,
                   (double) bench_result.perf_count / (n * bench_result.elt_per_op), (double) (n * bench_result.elt_per_op) / bench_result.perf_count,
                   bench_result.errors);
#           endif
        }
    }

    return 0;
}
