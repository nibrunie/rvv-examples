#include <math.h>
#include <stddef.h>


/** Baseline implementation of softmax for fp32
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
    for (i = 0; i < n; ++i) sum += expf(src[i]);

    // computing the reciprocal of the sum of exponentials, once and for all
    float inv_sum = 1.f / sum;

    // normalizing each element
    for (i = 0; i < n; ++i) dst[i] = src[i] / inv_sum;

}