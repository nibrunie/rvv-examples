#include <math.h>
#include <stddef.h>
#include <inttypes.h>
#include <riscv_vector.h>
#include <bench_softmax_utils.h>

softmax_bench_result_t softmax_bench(float* dst, float* src, softmax_func_t func, double* golden, size_t n); 


/** Quick and dirty implementation of exponential function on a binary32 input */
float quick_dirty_expf(float x) {
    const float ln2 = 0x1.62e43p-1;    
    const float iln2 = 0x1.715476p0f;

    // argument reduction
    const int k = nearbyintf(x * iln2);
    const float r = fmaf(- k, ln2, x);

    // polynomial approximation exp(r)
    // 0x1p0 + _x_ * (0x1.000002p0 +
    //         _x_ * (0x1.00001p-1 +
    //         _x_ * (0x1.55546ep-3 +
    //         _x_ * (0x1.554854p-5 +
    //         _x_ * (0x1.114662p-7 +
    //         _x_ * (0x1.7209d4p-10 +
    //         _x_ * 0x1.94480ap-13))))))
    const float poly_coeffs[] = {
        0x1.000002p0, 
        0x1.00001p-1, 
        0x1.55546ep-3, 
        0x1.554854p-5, 
        0x1.114662p-7, 
        0x1.7209d4p-10, 
        0x1.94480ap-13,
    };

    const int poly_degree = 6;

    float poly_r = poly_coeffs[poly_degree];
    int i = 0;
    for (i = poly_degree - 1; i >= 0; i--) {
        // poly_r = poly_r * r + poly_coeffs[i];
        poly_r = fmaf(poly_r, r, poly_coeffs[i]);
    }
    // poly_r = 1.f + r * poly_r;
    poly_r = fmaf(r, poly_r, 1.f);

    // NOTE: a proper cast should be done through memcopy and not an union
    typedef union { float f; uint32_t u; } f_u32_t;
    f_u32_t exp2_k;
    // quick and dirty (does not manage overflow/underflow/special values)
    // way to compute 2^k by injecting the biased exponent in the proper place
    // for IEEE-754 binary32 encoding.
    exp2_k.u = (127 + k) << 23;

    // result reconstruction
    float exp_x = poly_r * exp2_k.f;

    return exp_x;
}


/** Baseline implementation of softmax for binary32
 * 
 *  @param dst destination array
 *  @param src source array
 *  @param n   number of element(s)
*/
void softmax_rvv_fp32(float* dst, float* src, size_t n)
{
    int i;

    // computing the sum of exponentials
    float sum = 0.f;
    for (i = 0; i < n; ++i) sum += quick_dirty_expf(src[i]);//expf(src[i]);

    // computing the reciprocal of the sum of exponentials, once and for all
    float inv_sum = 1.f / sum;

    // normalizing each element
    size_t avl = n;
    while (avl > 0) {
        size_t vl = __riscv_vsetvl_e32m1(avl);
        vfloat32m1_t row = __riscv_vle32_v_f32m1(src, vl);
        row = __riscv_vfmul_vf_f32m1(row, inv_sum, vl);
        __riscv_vse32(dst, row, vl);
        avl -= vl;
        src += vl;
        dst += vl;
    }
}


softmax_bench_result_t softmax_rvv_fp32_bench(float* dst, float* src, double* golden, size_t n) {
    return softmax_bench(dst, src, softmax_rvv_fp32, golden, n);
}