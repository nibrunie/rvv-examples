#include <stddef.h>
#include <fenv.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <riscv_vector.h>

#include <bench_poly_mult_utils.h>

// function defined in ntt_scala.c
extern int ringPowers[1][8][64];
extern int ringInvPowers[1][8][64];

ring_t getRing(int degree);
void poly_fast_ntt_transform(ntt_t* dst, polynomial_t src, ring_t ring, int rootOfUnity);
void poly_fast_inv_ntt_tranform(polynomial_t* dst, ntt_t src, ring_t ring); 
void poly_fast_ntt_transform_helper(ntt_t* dst, int* coeffs, ring_t ring, int degree, int stride, int level, int rootPowers[8][64]);
void ntt_mul(ntt_t* dst, ntt_t lhs, ntt_t rhs); 

/** reinterpret cast help: reinterpreting to a different element size AND signedness */
static inline vint32m8_t __riscv_vreinterpret_v_u64m8_i32m8(vuint64m8_t vec) {
    vuint32m8_t vec_u32m8 = __riscv_vreinterpret_v_u64m8_u32m8(vec);
    return __riscv_vreinterpret_v_u32m8_i32m8(vec_u32m8);
}

/** Emulation of vwsll (when intrinsics is not defined) */
static inline vuint64m8_t __riscv_vwsll_vx_u64m8(vuint32m4_t vec, int amt, size_t vl) {
   return __riscv_vsll_vx_u64m8(__riscv_vzext_vf2_u64m8(vec, vl), amt, vl); 
}

void rvv_ntt_last_stage(ntt_t* dst, int* coeffs, int stride) {
    size_t avl = (dst->degree + 1) / 2;
    int* dst_coeffs = dst->coeffs;

    for (size_t vl; avl > 0; avl -= vl, coeffs += vl, dst_coeffs += 2*vl)
    {
        // compute loop body vector length from avl (application vector length)
        vl = __riscv_vsetvl_e32m4(avl);

        // loading even coefficients
        vuint32m4_t vec_even_coeffs = __riscv_vle32_v_u32m4((unsigned int*) coeffs, vl);
        // duplicating expansion
        vint32m8_t vec_lhs = __riscv_vreinterpret_v_u64m8_i32m8(
            __riscv_vor_vv_u64m8(
                __riscv_vzext_vf2_u64m8(vec_even_coeffs, vl),
                __riscv_vwsll_vx_u64m8(vec_even_coeffs, 32, vl),
                vl
            )
        );
        // loading odd coefficients
        vuint32m4_t vec_odd_coeffs = __riscv_vle32_v_u32m4((unsigned int*) coeffs + stride, vl);
        // duplicating expansion
        vint32m8_t vec_rhs = __riscv_vreinterpret_v_u64m8_i32m8(
            __riscv_vor_vv_u64m8(
                __riscv_vzext_vf2_u64m8(vec_odd_coeffs, vl),
                __riscv_vwsll_vx_u64m8(
                    __riscv_vreinterpret_v_i32m4_u32m4(
                        __riscv_vneg_v_i32m4(__riscv_vreinterpret_v_u32m4_i32m4(vec_odd_coeffs), vl)
                    ), 32, vl),
                vl
            )
        );
        vint32m8_t vec_rec = __riscv_vadd_vv_i32m8(
            vec_lhs,
            vec_rhs,
            vl*2); // FIXME: need cast to 32-bit

        // storing results
        __riscv_vse32_v_i32m8(dst_coeffs, vec_rec, 2*vl);
    }
}



/** Compute n-element NTT, assuming level @p level */
void rvv_ntt_stage(ntt_t* dst, int* coeffs, int n, int level, int rootPowers[8][64]) {
    const size_t coeffWidth = sizeof(coeffs[0]);

    size_t avl = n;
    ntt_t ntt_even = allocate_poly(n / 2, dst->modulo);
    ntt_t ntt_odd = allocate_poly(n / 2, dst->modulo);
    int* even_coeffs = ntt_even.coeffs;
    int* odd_coeffs = ntt_odd.coeffs; // coeffs + n / 2;
    for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl)
    {
        vl = __riscv_vsetvl_e32m8(avl);
        // splitting even and odd coefficients strided load
        vint32m8_t vec_even_coeffs = __riscv_vlse32_v_i32m8((int*) coeffs, sizeof(coeffs[0]), vl);
        vint32m8_t vec_odd_coeffs = __riscv_vlse32_v_i32m8((int*) (coeffs + 1), sizeof(coeffs[0]), vl);

        __riscv_vse32_v_i32m8((int*) even_coeffs, vec_even_coeffs, vl);
        __riscv_vse32_v_i32m8((int*) odd_coeffs, vec_odd_coeffs, vl);
    }

    // NTT recursion
    ntt_t dst_even = {.degree = (n / 2 - 1), .modulo = dst->modulo, coeffs = dst->coeffs, .coeffSize = (n/2)};
    ntt_t dst_odd = {.degree = (n / 2 - 1), .modulo = dst->modulo, coeffs = (dst->coeffs + n / 2), .coeffSize = (n/2)};
    rvv_ntt_stage(&dst_even, ntt_even.coeffs, n / 2, level + 1, rootPowers);
    rvv_ntt_stage(&dst_odd, ntt_odd.coeffs, n / 2, level + 1, rootPowers);

    // TODO: free ntt_even and ntt_odd

    // 
    avl = n;
    even_coeffs = dst->coeffs;
    odd_coeffs = dst->coeffs + n / 2;
    for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl)
    {
        vl = __riscv_vsetvl_e32m8(avl);
        vint32m8_t vec_even_coeffs = __riscv_vle32_v_i32m8((int*) even_coeffs, vl);
        vint32m8_t vec_odd_coeffs = __riscv_vle32_v_i32m8((int*) odd_coeffs, vl);

        vint32m8_t vec_twiddleFactor = __riscv_vle32_v_i32m8((int*) rootPowers[level], vl);

        // TODO: consider using a vectorized version of Barret's reduction algorithm
        vint32m8_t vec_odd_results = __riscv_vmul_vv_i32m8(vec_odd_coeffs, vec_twiddleFactor, vl);
        vint32m8_t vec_even_results = __riscv_vadd_vv_i32m8(vec_even_coeffs, vec_odd_results, vl);
        // even results
        vec_even_results = __riscv_vrem_vx_i32m8(vec_even_results, dst->modulo, vl);
        __riscv_vse32_v_i32m8(even_coeffs, vec_even_results, vl);

        // odd results
        vec_odd_results = __riscv_vsub_vv_i32m8(vec_even_coeffs, vec_odd_results, vl);
        vec_odd_results = __riscv_vrem_vx_i32m8(vec_odd_results, dst->modulo, vl);
        __riscv_vse32_v_i32m8(odd_coeffs, vec_odd_results, vl);
    }

}


void rvv_ntt_permute_inputs(ntt_t* dst, int* coeffs, int level) {
    // 2^level coefficients
    vint8m1_t mask_even_i8 = __riscv_vmv_v_x_i8m1(0x55, 16);
    vint8m1_t mask_odd_i8 = __riscv_vmv_v_x_i8m1(0xAA, 16);
    // mask type is bool<n> where n in EEW / EMUL (in our case 32 / 8 = 4)
    vbool4_t mask_even_b4 = __riscv_vreinterpret_v_i8m1_b4(mask_even_i8);
    vbool4_t mask_odd_b4 = __riscv_vreinterpret_v_i8m1_b4(mask_odd_i8);


    size_t avl = dst->degree + 1;
    size_t vl = __riscv_vsetvl_e32m8(avl);
    // TODO: loop around avl
    if (1) {
        // method 1: unit-strided load + vcompress
        // extracting even coefficients
        vuint32m8_t vec_even_coeffs = __riscv_vle32_v_u32m8((unsigned int*) coeffs, vl);
        vec_even_coeffs = __riscv_vcompress_vm_u32m8(vec_even_coeffs, mask_even_b4, vl);
        vuint32m8_t vec_odd_coeffs = __riscv_vle32_v_u32m8((unsigned int*) coeffs, vl);
        vec_odd_coeffs = __riscv_vcompress_vm_u32m8(vec_odd_coeffs, mask_odd_b4, vl);
    } else {
        // method 2: strided load
        vuint32m8_t vec_even_coeffs = __riscv_vlse32_v_u32m8((unsigned int*) coeffs, sizeof(coeffs[0]), vl);
        vuint32m8_t vec_odd_coeffs = __riscv_vlse32_v_u32m8((unsigned int*) (coeffs + 1), sizeof(coeffs[0]), vl);
    }
}


/**
 * @brief Multiply two NTT-transformed polynomials
 * @param dst destination polynomial
 * @param lhs left-hand-side polynomial
 * @param rhs right-hand-side polynomial
 *
 * @note The destination polynomial must have a degree greater or equal to the maximum of lhs and rhs degrees.
 */
void rvv_ntt_mul(ntt_t* dst, ntt_t lhs, ntt_t rhs) {
    assert(dst->degree >= lhs.degree && dst->degree >= rhs.degree);

    size_t avl = dst->degree + 1;
    int* lhs_coeffs = lhs.coeffs;
    int* rhs_coeffs = rhs.coeffs;
    int* dst_coeffs = dst->coeffs;
    for (size_t vl; avl > 0; avl -= vl, lhs_coeffs += vl, rhs_coeffs += vl, dst_coeffs += vl)
    {
        // compute loop body vector length from avl (application vector length)
        vl = __riscv_vsetvl_e32m8(avl);
        // loading operands
        vint32m8_t vec_src_lhs = __riscv_vle32_v_i32m8(lhs_coeffs, vl);
        vint32m8_t vec_src_rhs = __riscv_vle32_v_i32m8(rhs_coeffs, vl);
        // modulo multiplication (eventually we will want to consider other techniques
        // than a naive remainder; e.g. Barret's reduction algorithm using a pre-computed
        // factor from the static modulo).
        vint32m8_t vec_acc = __riscv_vmul_vv_i32m8(vec_src_lhs, vec_src_rhs, vl);
        // TODO: consider using a vectorized version of Barret's reduction algorithm
        vec_acc = __riscv_vrem_vx_i32m8(vec_acc, dst->modulo, vl);
        // storing results
        __riscv_vse32_v_i32m8(dst_coeffs, vec_acc, vl);
    }
}

/** dividing each coefficient by the degree */
void rvv_ntt_degree_scaling(ntt_t* dst, ring_t ring) {
    size_t avl = dst->degree + 1;
    int* dst_coeffs = dst->coeffs;
    for (size_t vl; avl > 0; avl -= vl, dst_coeffs += vl)
    {
        // compute loop body vector length from avl (application vector length)
        vl = __riscv_vsetvl_e32m8(avl);
        // loading operands
        vint32m8_t vec = __riscv_vle32_v_i32m8(dst_coeffs, vl);
        // modulo multiplication (eventually we will want to consider other techniques
        // than a naive remainder; e.g. Barret's reduction algorithm using a pre-computed
        // factor from the static modulo).
        vec = __riscv_vmul_vx_i32m8(vec, ring.invDegree, vl);
        // modulo reduction
        vec = __riscv_vrem_vx_i32m8(vec, ring.modulo, vl);
        // storing results
        __riscv_vse32_v_i32m8(dst_coeffs, vec, vl);
    }
}

void poly_mult_ntt_rvv(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(dst->degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    poly_fast_ntt_transform_helper(&ntt_lhs, lhs.coeffs, ring, lhs.degree, 1, 0, ringPowers[0]);
    poly_fast_ntt_transform_helper(&ntt_lhs_times_rhs, rhs.coeffs, ring, rhs.degree, 1, 0, ringPowers[0]);

    // element-size multiplication using RVV
    rvv_ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    poly_fast_ntt_transform_helper(dst, ntt_lhs_times_rhs.coeffs, ring, ntt_lhs_times_rhs.degree, 1, 0, ringInvPowers[0]);
    // division by the degree
    rvv_ntt_degree_scaling(dst, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

/** Benchmark wrapper */
poly_mult_bench_result_t poly_mult_ntt_rvv_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv, golden);
}