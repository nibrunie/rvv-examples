#include <math.h>
#include <stddef.h>
#include <bench_softmax_utils.h>


/** Baseline implementation of softmax for binary32
 * 
 *  @param dst destination array
 *  @param src source array
 *  @param n   number of element(s)
*/
void softmax_baseline_fp32(float* dst, float* src, size_t n)
{
    int i;

    // computing the sum of exponentials
    float sum = 0.f;
    for (i = 0; i < n; ++i) {
        dst[i] = expf(src[i]);
        sum += dst[i]; 
    }

    // computing the reciprocal of the sum of exponentials, once and for all
    float inv_sum = 1.f / sum;

    // normalizing each element
    for (i = 0; i < n; ++i) dst[i] = dst[i] * inv_sum;
}


float quick_dirty_expf(float x); 

void softmax_scalar_quick_dirty_expf_fp32(float* dst, float* src, size_t n)
{
    int i;

    // computing the sum of exponentials
    float sum = 0.f;
    for (i = 0; i < n; ++i) {
        dst[i] = quick_dirty_expf(src[i]);
        sum += dst[i]; 
    }

    // computing the reciprocal of the sum of exponentials, once and for all
    float inv_sum = 1.f / sum;

    // normalizing each element
    for (i = 0; i < n; ++i) dst[i] = dst[i] * inv_sum;
}

/** Baseline implementation of softmax on binary32 input using binary64
 * 
 *  @param dst destination array
 *  @param src source array
 *  @param n   number of element(s)
*/
void softmax_golden_fp32_fp64(double* dst, float* src, size_t n)
{
    int i;
    // looking for max input value
    double max_x = src[0];
    for (i = 1; i < n; ++i) {
        if (src[i] > max_x) max_x = src[i]; 
    }

    // computing the sum of exponentials
    double sum = 0.;
    for (i = 0; i < n; ++i) {
        dst[i] = exp((double) src[i] - max_x);
        sum += dst[i];
    }

    // normalizing each element: we use the division every single
    // type rather than hoisting the evaluation of the reciprocal 1 / sum
    // to improve accuracy since this function is intended to build golden
    // values and not to be fast.
    for (i = 0; i < n; ++i) dst[i] = (double) dst[i] / sum;
}

/** generic benchmark wrapper for n-element 1D softmax implementation
 *
 * @param dst destination array (must be at least n-element wide)
 * @param src source array (must be at least n-element wide)
 * @param func function under test
 * @param goldem golden value for destination array (must contain n valid elements)
 * @param n number of elements to be considered
*/
softmax_bench_result_t softmax_bench(float* dst, float* src, softmax_func_t func, double* golden, size_t n) {
    unsigned long start, stop;
    softmax_bench_result_t bench_result;

    start = read_perf_counter();
    func(dst, src, n);
    stop = read_perf_counter();
    bench_result.perf_count = stop - start;

    // initializing bench result error values
    bench_result.max_abs_error = 0.0;
    bench_result.max_rel_error = 0.0;

    unsigned i;
    double sum_square_errors = 0.0;
    double sum_rel_errors = 0.0;
#       ifdef VERY_VERBOSE
        printf("relative errors: ");
#       endif
    for (i = 0; i < n; ++i) {
        double abs_error = fabs(dst[i] - golden[i]);
        double rel_error = abs_error / (double) golden[i];
#       ifdef VERY_VERBOSE
        printf("%.8e ", rel_error);
#       endif
        sum_square_errors += rel_error * rel_error;
        sum_rel_errors += rel_error;

        if (abs_error > bench_result.max_abs_error) bench_result.max_abs_error = abs_error;
        if (rel_error > bench_result.max_rel_error) bench_result.max_rel_error = rel_error;
    }

#   ifdef VERY_VERBOSE
    printf("\n");
#   endif

    bench_result.error_norm2 = sqrt(sum_square_errors);
    bench_result.mean_rel_error = sum_rel_errors / (double) n;
    return bench_result;
}


softmax_bench_result_t softmax_baseline_fp32_bench(float* dst, float* src, double* golden, size_t n) {
    return softmax_bench(dst, src, softmax_baseline_fp32, golden, n);
}


softmax_bench_result_t softmax_scalar_quick_dirty_expf_fp32_bench(float* dst, float* src, double* golden, size_t n) {
    return softmax_bench(dst, src, softmax_scalar_quick_dirty_expf_fp32, golden, n);
}