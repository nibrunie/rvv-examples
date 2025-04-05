#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "bench_utils.h"

// setting a default message size (1MiB) if none is defined
#ifndef MSG_SIZES
# define MSG_SIZES (1024 * 1024)
#endif 


int32_t golden_min(int32_t* vec, size_t len) {
    if (len <= 0) return -1;
    int32_t min = vec[0];
    for (size_t i = 1; i < len; ++i) {
        if (vec[i] < min) {
            min = vec[i];
        }
    }
    return min;
}

int32_t rvv_min(int32_t* vec, size_t len) {
    int32_t min = INT32_MAX;
    asm __volatile__ (
        "vsetivli x0, 1, e32, m1\n"     // set vector length, before initializing the scalar value
        "vmv.v.x v16, %[min]\n"         // initializing the minimum value in v8 (vector register)
    "rvv_min_loop:\n"                              // label for the loop
        "vsetvli a4, %[avl], e32, m8\n" // set vector length
        "vle32.v v8, (%[vec])\n"        // load vector from memory
        "vredmin.vs v16, v8, v16\n"     // compute the minimum in the vector
        "sub %[avl], %[avl], a4\n"      // decrement the remaining length
        "sh2add %[vec], a4, %[vec]\n"   // move the vector pointer forward by the number of processed elements
        "bnez %[avl], rvv_min_loop\n"             // if there are more elements to process, loop back
        "vmv.x.s %[min], v16\n"         // move the final minimum value from vector register to scalar
        : [avl]"+r"(len), [min]"+r"(min), [vec]"+r"(vec)                 // output operand (min)
        :
        : "memory", "a4", "v8", "v16"
    );
    return min;
}


int32_t rvv_min_2_steps(int32_t* vec, size_t len) {
    int32_t min = INT32_MAX;
    asm __volatile__ (
        "vsetvli a4,x0, e32, m8\n"    // set vector length, before initializing the scalar value
        "vmv.v.x v16, %[min]\n"         // initializing the minimum value in v8 (vector register)
    "1:\n"                   // label for the loop
        "vsetvli a4, %[avl], e32, m8\n" // set vector length
        "vle32.v v8, (%[vec])\n"        // load vector from memory
        "vmin.vv v16, v8, v16\n"        // compute the minimum in the vector
        "sub %[avl], %[avl], a4\n"      // decrement the remaining length
        "sh2add %[vec], a4, %[vec]\n"   // move the vector pointer forward by the number of processed elements
        "bnez %[avl], 1b\n"   // if there are more elements to process, loop back

        "vsetvli a4, x0, e32, m8\n"    // set vector length, before initializing the scalar value
        "vredmin.vs v16, v16, v16\n"    // this is a second step to ensure the final reduction in case of multiple vector registers
        "vmv.x.s %[min], v16\n"         // move the final minimum value from vector register to scalar
        : [avl]"+r"(len), [min]"+r"(min), [vec]"+r"(vec)                 // output operand (min)
        :
        : "memory", "a4", "v8", "v16"
    );
    return min;
}

float vector_dot_product(float* lhs, float* rhs, size_t len) {
    float result = 0.0f;
    for (size_t i = 0; i < len; ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}


typedef struct {
    int32_t (*red_func)(int32_t* vec, size_t len);
    int32_t (*golden)(int32_t* vec, size_t len);
    char* label;
} reduction_bench_t;

int bench_int_reduction(void) {
    uint32_t start = 0, stop = 0;

    reduction_bench_t benchmarks[] = {
        // labels should be padded to all have the same width (result display alignment)
        {.red_func = golden_min,      .golden = golden_min, .label = "generic vector min reduction"},
        {.red_func = rvv_min,         .golden = golden_min, .label = "RVV based vector min reduction"},
        {.red_func = rvv_min_2_steps, .golden = golden_min, .label = "RVV based 2-step vector min reduction"},
    };

    size_t vectorSizes[] = {16, 32, 64, 128, 512, 1024, 16384, MSG_SIZES};
    int error = 0;
#   ifndef VERBOSE
        printf("message_size, red_result, label, perf(" PERF_METRIC ")\n");
#   endif

    for (int size_id = 0; size_id < sizeof(vectorSizes) / sizeof(size_t); ++size_id) {
        size_t vecSize = vectorSizes[size_id];
        int32_t* inputVec = (int32_t*) malloc(vecSize * sizeof(int32_t)); // allocate space for the vector
        assert(inputVec != NULL && "Failed to allocate memory for input vector"); // ensure allocation was successful

        for (size_t i = 0; i < vecSize; ++i) {
            inputVec[i] = (int32_t)(rand()); // fill the vector with random integers in the range [0, 100)
        }

        for (int bench_id = 0; bench_id < sizeof(benchmarks) / sizeof(reduction_bench_t); bench_id++) {
            start = read_perf_counter();
            uint32_t result = benchmarks[bench_id].red_func(inputVec, vecSize); // call the reduction function
            stop = read_perf_counter();

            // computing golden
            int32_t golden = benchmarks[bench_id].golden(inputVec, vecSize);

            // checks
            error += result != golden;

            uint32_t delay = stop - start;
            float throughput = (double) delay / vecSize; 

            // full message
    #       ifdef VERBOSE
            printf("reduction(vector[%lu]) = %"PRIi32" (%s)      in %u " PERF_METRIC "(s) [%.3f " PERF_METRIC "(s) per Byte]\n",
                vecSize, result, benchmarks[bench_id].label, delay, throughput);
    #       else
            printf("%8lu, %8"PRIi32", %s, %u\n", vecSize, result, benchmarks[bench_id].label, delay);
    #       endif

        }

        free(inputVec);

    }

    return error;
}


int main(void) {
    return bench_int_reduction();
}