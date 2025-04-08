#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "bench_utils.h"

// setting a default message size (1MiB) if none is defined
#ifndef VECTOR_SIZE
# define VECTOR_SIZE (1024 * 1024)
#endif 

#define xstr(x) #x
#define str(x) xstr(x)

#ifndef LMUL
#define LMUL m8
#endif

/** Golden reference implementation of vector reduction-based minimum
 *  of 32-bit signed integers.
 *
 * @param vec input vector
 * @param len number of elements in the input vector
 * @return minimum of the vector element(s)
 */
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

#define BENCH_REDUCTION_OP(op, init_value) \
__attribute__((noinline)) int32_t reduction_insn_bench_ ## op (int32_t* vec, size_t len) { \
    int32_t acc = init_value; \
    asm __volatile__ ( \
        "vsetivli x0, 1, e32, m1, ta, ma \n" \
        "vmv.v.x v16, %[acc]\n" \
        "vsetvli a4, %[avl], e32, " str(LMUL) ", tu, mu\n" \
	"vle32.v v8, (%[vec])\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        str(op) ".vs v16, v8, v16\n" \
        "vmv.x.s %[acc], v16\n" \
        : [avl]"+r"(len), [acc]"+r"(acc), [vec]"+r"(vec) \
        : \
        : "a4", "v8", "v16" \
    ); \
    return acc; \
}

BENCH_REDUCTION_OP(vredmin, INT32_MAX) // Generate the reduction instruction for minimum
BENCH_REDUCTION_OP(vredsum, 0) // Generate the reduction instruction for sum (for completeness)
BENCH_REDUCTION_OP(vfredusum, 0.f)
BENCH_REDUCTION_OP(vfredosum, 0.f) // Generate the reduction instruction for ordered sum
				       //

typedef struct {
    int32_t (*red_func)(int32_t* vec, size_t len);
    char* label;
} synthetic_reduction_bench_t;

int synthetic_bench() {
    int vlmax = -1;
    __asm __volatile__ (
        "vsetvli %[vlmax], x0, e32, " str(LMUL) ", ta, ma \n" // set vector length
        : [vlmax]"=r"(vlmax)
        :
        : "memory"
    );

    synthetic_reduction_bench_t benchmarks[] = {
        // labels should be padded to all have the same width (result display alignment)
        {.red_func = reduction_insn_bench_vredmin, .label = "vredmin.vs;lmul=" str(LMUL)}, // this is the minimum reduction
        {.red_func = reduction_insn_bench_vredsum, .label = "vredsum.vs;lmul=" str(LMUL)}, // sum reduction
        {.red_func = reduction_insn_bench_vfredusum, .label = "vfredusum.vs;lmul=" str(LMUL)}, // unordered sum
        {.red_func = reduction_insn_bench_vfredosum, .label = "vfredosum.vs;lmul=" str(LMUL)}, // ordered sum
    };

#   ifndef VERBOSE
        printf("vector_size, result, label, perf(" PERF_METRIC ")\n");
#   endif

    int32_t* inputVec = (int32_t*) malloc(vlmax * sizeof(int32_t)); // allocate space for the vector
    for (int i = 0; i <= vlmax; i++) inputVec[i] = 0x3f800017 + (i << 8) + (i % 3) + ((1u % 2) << 31) + ((i % 5) << 25);
    uint32_t start = 0, stop = 0;

    for (int vl = 0; vl <= vlmax; ++vl) {
        for (int bench_id = 0; bench_id < sizeof(benchmarks) / sizeof(synthetic_reduction_bench_t); bench_id++) {
            start = read_perf_counter();
            int32_t result = benchmarks[bench_id].red_func(inputVec, vl); // call the reduction function
            stop = read_perf_counter();

	    // dividing by 30 because we run each instruction 30 time
            float delay = (stop - start) / 30.f;
            float throughput = (double) delay / vl; 

            // full message
    #       ifdef VERBOSE
            printf("reduction(vector[%lu]) = %"PRIi32" (%s)      in %.2f " PERF_METRIC "(s) [%.3f " PERF_METRIC "(s) per Byte]\n",
                vl, result, benchmarks[bench_id].label, delay, throughput);
    #       else
            printf("%8lu, %8"PRIx32", %s, %.2f\n", vl, result, benchmarks[bench_id].label, delay);
    #       endif
        }
    }

    free(inputVec); // free the allocated vector
    return 0; // return success

}


/** RVV-based implementation of vector minimum
 *  of 32-bit signed integers.
 *
 * @param vec input vector
 * @param len number of elements in the input vector
 * @return minimum of the vector element(s)
 */
__attribute__((noinline)) int32_t rvv_min(int32_t* vec, size_t len) {
    int32_t min = INT32_MAX;
    asm __volatile__ (
        "vsetivli x0, 1, e32, m1, ta, ma \n"     // set vector length, before initializing the scalar value
        "vmv.v.x v16, %[min]\n"         // initializing the minimum value in v8 (vector register)
    "rvv_min_loop:\n"                              // label for the loop
        "vsetvli a4, %[avl], e32, " str(LMUL) ", ta, ma\n" // set vector length
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
        "vsetvli a4,x0, e32, " str(LMUL) ", ta, ma\n"    // set vector length, before initializing the scalar value
        "vmv.v.x v16, %[min]\n"         // initializing the minimum value in v8 (vector register)
    "1:\n"                   // label for the loop
        "vsetvli a4, %[avl], e32, " str(LMUL) ", ta, ma\n" // set vector length
        "vle32.v v8, (%[vec])\n"        // load vector from memory
        "vmin.vv v16, v8, v16\n"        // compute the minimum in the vector
        "sub %[avl], %[avl], a4\n"      // decrement the remaining length
        "sh2add %[vec], a4, %[vec]\n"   // move the vector pointer forward by the number of processed elements
        "bnez %[avl], 1b\n"   // if there are more elements to process, loop back

        "vsetvli a4, x0, e32, " str(LMUL) ", ta, ma\n"    // set vector length, before initializing the scalar value
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


/** RVV-based vector dot product using ordered reduction sum
 *
 * @param lhs pointer to the left hand side vector
 * @param rhs pointer to the right hand side vector
 * @param len length of the vectors
 * @return the dot product of the two vectors
 */
__attribute__((noinline)) float rvv_dot_product(float* lhs, float* rhs, size_t len) {
    float res = 0.f;
    asm __volatile__ (
        "vsetivli x0, 1, e32, m1, ta, ma\n"     // set vector length, before initializing the scalar value
        "vfmv.v.f v24, %[res]\n"         // initializing the temporary accumulator in v24[0] to zero
    "rvv_dot_product_loop:\n"                              // label for the loop
        "vsetvli a4, %[avl], e32, " str(LMUL) ", ta, ma\n" // set vector length
        "vle32.v v8, (%[lhs])\n"        // load left hand side vector from memory
        "vle32.v v16, (%[rhs])\n"        // load right hand side vector from memory
        "vfmul.vv v8, v8, v16\n"         // multiply the two vectors element-wise
        "vfredosum.vs v24, v8, v24\n"     // accumulate the product into the temporary accumulator (v24[0])
        "sub %[avl], %[avl], a4\n"      // decrement the remaining length
        "sh2add %[lhs], a4, %[lhs]\n"   // move the vector pointer forward by the number of processed elements
        "sh2add %[rhs], a4, %[rhs]\n"   // move the vector pointer forward by the number of processed elements
        "bnez %[avl], rvv_dot_product_loop\n"             // if there are more elements to process, loop back
        "vfmv.f.s %[res], v24\n"         // move the final minimum value from vector register to scalar
        : [avl]"+r"(len), [res]"+f"(res), [lhs]"+r"(lhs), [rhs]"+r"(rhs)
        :
        : "memory", "a4", "v8", "v16", "v24" // clobbered registers
    );
    return res;
}


/** RVV-based vector dot product using unordered reduction sum
 *
 * @param lhs pointer to the left hand side vector
 * @param rhs pointer to the right hand side vector
 * @param len length of the vectors
 * @return the dot product of the two vectors
 */
__attribute__((noinline)) float rvv_dot_product_u(float* lhs, float* rhs, size_t len) {
    float res = 0.f;
    asm __volatile__ (
        "vsetivli x0, 1, e32, m1, ta, ma\n"     // set vector length, before initializing the scalar value
        "vfmv.v.f v24, %[res]\n"         // initializing the temporary accumulator in v24[0] to zero
    "rvv_dot_product_u_loop:\n"                              // label for the loop
        "vsetvli a4, %[avl], e32, " str(LMUL) ", ta, ma\n" // set vector length
        "vle32.v v8, (%[lhs])\n"        // load left hand side vector from memory
        "vle32.v v16, (%[rhs])\n"        // load right hand side vector from memory
        "vfmul.vv v8, v8, v16\n"         // multiply the two vectors element-wise
        "vfredusum.vs v24, v8, v24\n"     // accumulate the product into the temporary accumulator (v24[0])
        "sub %[avl], %[avl], a4\n"      // decrement the remaining length
        "sh2add %[lhs], a4, %[lhs]\n"   // move the vector pointer forward by the number of processed elements
        "sh2add %[rhs], a4, %[rhs]\n"   // move the vector pointer forward by the number of processed elements
        "bnez %[avl], rvv_dot_product_u_loop\n"             // if there are more elements to process, loop back
        "vfmv.f.s %[res], v24\n"         // move the final minimum value from vector register to scalar
        : [avl]"+r"(len), [res]"+f"(res), [lhs]"+r"(lhs), [rhs]"+r"(rhs)
        :
        : "memory", "a4", "v8", "v16", "v24" // clobbered registers
    );
    return res;
}

/** RVV-based vector dot product using element wise accumulation
 *
 * @param lhs pointer to the left hand side vector
 * @param rhs pointer to the right hand side vector
 * @param len length of the vectors
 * @return the dot product of the two vectors
 */
__attribute__((noinline)) float rvv_dot_product_2_steps(float* lhs, float* rhs, size_t len) {
    float res = 0.f;
    asm __volatile__ (
        "vsetvli a4, x0, e32, " str(LMUL) ", ta, ma\n" // set vector length
        "vfmv.v.f v24, %[res]\n"    // initializing the temporary accumulators in v24 to zeros
        "vmv.v.i v24, 0\n"
    "rvv_dot_product_2_steps_loop:\n"                   // label for the loop
        "vsetvli a4, %[avl], e32, " str(LMUL) ", tu, mu\n"         // set vector length
        "vle32.v v8, (%[lhs])\n"                        // load left hand side vector from memory
        "vle32.v v16, (%[rhs])\n"                       // load right hand side vector from memory
        "vfmacc.vv v24, v8, v16\n"                      // multiply the two vectors element-wise
        "sub %[avl], %[avl], a4\n"                      // decrement the remaining length
        "sh2add %[lhs], a4, %[lhs]\n"                   // move the vector pointer forward by the number of processed elements
        "sh2add %[rhs], a4, %[rhs]\n"                   // move the vector pointer forward by the number of processed elements
        "bnez %[avl], rvv_dot_product_2_steps_loop\n"   // if there are more elements to process, loop back

        "vsetvli a4, x0, e32, " str(LMUL) ", ta, ma\n" // set vector length
        "vmv.v.i v16, 0\n" // initialize v16 to zero for the final reduction step
        "vfredusum.vs v16, v24, v16\n"   // this is a second step to ensure the final reduction in case of multiple vector registers
        "vfmv.f.s %[res], v16\n"         // move the final minimum value from vector register to scalar
        : [avl]"+r"(len), [res]"+f"(res), [lhs]"+r"(lhs), [rhs]"+r"(rhs)
        :
        : "memory", "a4", "v8", "v16", "v24" // clobbered registers
    );
    return res;
}


/** Integer reduction benchmark (single array interface)
 *
 * @param red_func the reduction function to benchmark
 * @param golden the golden reference implementation for comparison
 * @param label a label describing the benchmark
 */
typedef struct {
    int32_t (*red_func)(int32_t* vec, size_t len);
    int32_t (*golden)(int32_t* vec, size_t len);
    char* label;
} int_reduction_bench_t;

/** Floating-point reduction benchmark (2-array interface)
 *
 * @param red_func the reduction function to benchmark
 * @param golden the golden reference implementation for comparison
 * @param label a label describing the benchmark
 */
typedef struct {
    float (*red_func)(float* lhs, float*rhs, size_t len);
    float (*golden)(float* lhs, float* rhs, size_t len);
    char* label;
} float_reduction_bench_t;


/** Integer reduction benchmark
 *
 * @return zero if all benchmarks passed, otherwise returns the number of failed benchmarks.
 */
int bench_int_reduction(void) {
    uint32_t start = 0, stop = 0;

    int_reduction_bench_t benchmarks[] = {
        // labels should be padded to all have the same width (result display alignment)
        {.red_func = golden_min,      .golden = golden_min, .label = "generic vector min reduction"},
        {.red_func = rvv_min,         .golden = golden_min, .label = "RVV based vector min reduction"},
        {.red_func = rvv_min_2_steps, .golden = golden_min, .label = "RVV based 2-step vector min reduction"},
    };

    size_t vectorSizes[] = {16, 32, 64, 128, 512, 1024, 16384, VECTOR_SIZE};
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

        for (int bench_id = 0; bench_id < sizeof(benchmarks) / sizeof(int_reduction_bench_t); bench_id++) {
            start = read_perf_counter();
            int32_t result = benchmarks[bench_id].red_func(inputVec, vecSize); // call the reduction function
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


/** Floating-point reduction benchmark
 *
 * @return zero if all benchmarks passed, otherwise returns the number of failed benchmarks.
 */
int bench_float_reduction(void) {
    uint32_t start = 0, stop = 0;

    float_reduction_bench_t benchmarks[] = {
        // labels should be padded to all have the same width (result display alignment)
        {.red_func = vector_dot_product,      .golden = vector_dot_product, .label = "generic vector dot product"},
        {.red_func = rvv_dot_product,         .golden = vector_dot_product, .label = "RVV-based vector dot product (ordered vfredosum)"},
        {.red_func = rvv_dot_product_u,       .golden = vector_dot_product, .label = "RVV-based vector dot product (unordered vfredusum)"},
        {.red_func = rvv_dot_product_2_steps, .golden = rvv_dot_product, .label = "RVV-based vector dot product (parallel accumulation)"},
    };

    size_t vectorSizes[] = {16, 32, 64, 128, 511, 512, 513, 1024, 16384, VECTOR_SIZE};
    int error = 0;
#   ifndef VERBOSE
        printf("message_size, red_result, label, perf(" PERF_METRIC ")\n");
#   endif

    for (int size_id = 0; size_id < sizeof(vectorSizes) / sizeof(size_t); ++size_id) {
        size_t vecSize = vectorSizes[size_id];
        float* inputVecLHS = (float*) malloc(vecSize * sizeof(float));
        assert(inputVecLHS != NULL && "Failed to allocate memory for input vector"); // ensure allocation was successful
        float* inputVecRHS = (float*) malloc(vecSize * sizeof(float));
        assert(inputVecRHS != NULL && "Failed to allocate memory for input vector"); // ensure allocation was successful

        for (size_t i = 0; i < vecSize; ++i) {
            // fill the input vectors with random floats in the range [0, 7)
            // floats are actually converted integer so we do not get
            // rounding errors for small enough vector sizes.
            inputVecLHS[i] = (float)(i % 7);
            inputVecRHS[i] = (float)(i % 7);
        }

        for (int bench_id = 0; bench_id < sizeof(benchmarks) / sizeof(float_reduction_bench_t); bench_id++) {
            start = read_perf_counter();
            float result = benchmarks[bench_id].red_func(inputVecLHS, inputVecRHS, vecSize); // call the reduction function
            stop = read_perf_counter();

            // computing golden
            float golden = benchmarks[bench_id].golden(inputVecLHS, inputVecRHS, vecSize);

            // checks
            error += result != golden;

            uint32_t delay = stop - start;
            float throughput = (double) delay / vecSize; 

            // full message
    #       ifdef VERBOSE
            printf("reduction(vector[%lu]) = %f (%s)      in %u " PERF_METRIC "(s) [%.3f " PERF_METRIC "(s) per Byte]\n",
                vecSize, result, benchmarks[bench_id].label, delay, throughput);
    #       else
            printf("%8lu, %f, %s, %u\n", vecSize, result, benchmarks[bench_id].label, delay);
    #       endif

        }
        free(inputVecLHS);
        free(inputVecRHS);
    }

    return error;
}


int main(void) {
    int error = 0;
    
    error |= bench_int_reduction();

    error |= bench_float_reduction();

    error |= synthetic_bench();

    return error;
}
