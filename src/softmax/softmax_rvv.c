#include <math.h>
#include <stddef.h>
#include <inttypes.h>
#include <riscv_vector.h>
#include <bench_softmax_utils.h>

softmax_bench_result_t softmax_bench(float* dst, float* src, softmax_func_t func, double* golden, size_t n); 


/** Quick and dirty implementation of exponential function on a binary32 input */
float quick_dirty_expf(float x) {
    // values determined using (python)sollya
    // >>> iln2 = sollya.round(1/sollya.log(2), sollya.binary32, sollya.RN)
    // >>> ln2 = sollya.round(sollya.log(2), sollya.binary32, sollya.RN)
    const float ln2 = 0x1.62e43p-1;    
    const float iln2 = 0x1.715476p0f;

    // argument reduction
    const int k = nearbyintf(x * iln2);
    const float r = fmaf(- k, ln2, x);

    // polynomial approximation exp(r)
    // coefficients determined using (python)sollya
    // >>> ln2ov2 = sollya.round(sollya.log(2), sollya.binary32, sollya.RN)
    // >>> approxInt = sollya.Interval(-ln2ov2, ln2ov2)
    // >>> approxFun = sollya.exp(sollya.x)
    // >>> degree = 7
    // >>> poly = sollya.fpminimax(approxFunc,
    //                             degree,
    //                             [1] + [sollya.binary32] * degree,
    //                             approxInterval)
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

/** RVV-based vectorized implementation of binary32 exponential with 
 *  result reduction (sum).
*/
float quick_dirty_vector_expf(float* dst, float* src, size_t n) {
    // values determined using (python)sollya
    const float ln2 = 0x1.62e43p-1;    
    const float iln2 = 0x1.715476p0f;

    const size_t vlmax = __riscv_vsetvlmax_e32m1(); 
    const vfloat32m1_t vln2 = __riscv_vfmv_v_f_f32m1(ln2, vlmax);
    const vfloat32m1_t viln2 = __riscv_vfmv_v_f_f32m1(iln2, vlmax);

    // element-wise reduction accumulator
    vfloat32m1_t vsum = __riscv_vfmv_v_f_f32m1(0.f, vlmax);

    size_t avl = n;
    while (avl > 0) {
        size_t vl = __riscv_vsetvl_e32m1(avl);
        vfloat32m1_t vx = __riscv_vle32_v_f32m1(src, vl);

        // argument reduction
        vfloat32m1_t vxiln2 = __riscv_vfmul(vx, iln2, vl);
        vint32m1_t       vk = __riscv_vfcvt_x_f_v_i32m1(vxiln2, vl); // require round to nearest mode
        vfloat32m1_t    vfk = __riscv_vfcvt_f_x_v_f32m1(vk, vl);
        // using vfnmsac.vf to evaluate r = x - k * log(2)
        vfloat32m1_t     vr = __riscv_vfnmsac(vx, ln2, vfk, vl);

        // polynomial approximation exp(r)
        const float poly_coeffs[] = {
            0x1.000002p0, 
            0x1.00001p-1, 
            0x1.55546ep-3, 
            0x1.554854p-5, 
            0x1.114662p-7, 
            0x1.7209d4p-10, 
            0x1.94480ap-13,
        };

        vfloat32m1_t poly_vr = __riscv_vfmv_v_f_f32m1(poly_coeffs[6], vl);
        poly_vr = __riscv_vfmadd(poly_vr, vr, __riscv_vfmv_v_f_f32m1(poly_coeffs[5], vl), vl);
        poly_vr = __riscv_vfmadd(poly_vr, vr, __riscv_vfmv_v_f_f32m1(poly_coeffs[4], vl), vl);
        poly_vr = __riscv_vfmadd(poly_vr, vr, __riscv_vfmv_v_f_f32m1(poly_coeffs[3], vl), vl);
        poly_vr = __riscv_vfmadd(poly_vr, vr, __riscv_vfmv_v_f_f32m1(poly_coeffs[2], vl), vl);
        poly_vr = __riscv_vfmadd(poly_vr, vr, __riscv_vfmv_v_f_f32m1(poly_coeffs[1], vl), vl);
        poly_vr = __riscv_vfmadd(poly_vr, vr, __riscv_vfmv_v_f_f32m1(poly_coeffs[0], vl), vl);
        poly_vr = __riscv_vfmadd(poly_vr, vr, __riscv_vfmv_v_f_f32m1(1.f, vl), vl);

        // reconstruction
        const int exp_bias = 127;
        vint32m1_t vbiased_exp = __riscv_vadd(vk, exp_bias, vl);
        vint32m1_t vexp2_vk    = __riscv_vsll(vbiased_exp, 23, vl);
        vfloat32m1_t vfexp2_vk = __riscv_vreinterpret_v_i32m1_f32m1(vexp2_vk);

        vfloat32m1_t vexp_vx  = __riscv_vfmul(poly_vr, vfexp2_vk, vl);

        // element-size reduction with redution accumulator
        // tail-undisturbed is mandatory here to ensure that if vl is less
        // than VLMAX then unaffacted sum terms are not changed.
        vsum = __riscv_vfadd_vv_f32m1_tu(vsum, vsum, vexp_vx, vl);

        __riscv_vse32(dst, vexp_vx, vl);
        avl -= vl;
        src += vl;
        dst += vl;
    }

    vfloat32m1_t vredsum = __riscv_vfmv_v_f_f32m1(0.f, vlmax);
    vredsum = __riscv_vfredosum_vs_f32m1_f32m1(vsum, vredsum, vlmax);

    return __riscv_vfmv_f_s_f32m1_f32(vredsum);
}


/** implementation of softmax for binary32 with RVV-based normalization
 * 
 *  @param dst destination array
 *  @param src source array
 *  @param n   number of element(s)
*/
void softmax_rvv_norm_fp32(float* dst, float* src, size_t n)
{
    int i;

    // computing the sum of exponentials
    float sum = 0.f;
    for (i = 0; i < n; ++i) {
        dst[i] = quick_dirty_expf(src[i]); //expf(src[i]);
        sum += dst[i];
    }

    // computing the reciprocal of the sum of exponentials, once and for all
    float inv_sum = 1.f / sum;

    // normalizing each element
    size_t avl = n;
    while (avl > 0) {
        size_t vl = __riscv_vsetvl_e32m1(avl);
        vfloat32m1_t row = __riscv_vle32_v_f32m1(dst, vl);
        row = __riscv_vfmul(row, inv_sum, vl);
        __riscv_vse32(dst, row, vl);
        avl -= vl;
        dst += vl;
    }
}


softmax_bench_result_t softmax_rvv_norm_fp32_bench(float* dst, float* src, double* golden, size_t n) {
    return softmax_bench(dst, src, softmax_rvv_norm_fp32, golden, n);
}


/** implementation of softmax for binary32 RVV-based
 * 
 *  @param dst destination array
 *  @param src source array
 *  @param n   number of element(s)
*/
void softmax_rvv_fp32(float* dst, float* src, size_t n)
{
    // computing element-wise exponentials and their seum
    float sum = quick_dirty_vector_expf(dst, src, n);

    // computing the reciprocal of the sum of exponentials, once and for all
    float inv_sum = 1.f / sum;

    // normalizing each element
    size_t avl = n;
    while (avl > 0) {
        size_t vl = __riscv_vsetvl_e32m1(avl);
        vfloat32m1_t row = __riscv_vle32_v_f32m1(dst, vl);
        row = __riscv_vfmul_vf_f32m1(row, inv_sum, vl);
        __riscv_vse32(dst, row, vl);
        avl -= vl;
        dst += vl;
    }
}


softmax_bench_result_t softmax_rvv_fp32_bench(float* dst, float* src, double* golden, size_t n) {
    return softmax_bench(dst, src, softmax_rvv_fp32, golden, n);
}