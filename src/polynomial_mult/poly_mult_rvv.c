#include <stddef.h>
#include <fenv.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <riscv_vector.h>

#include <bench_poly_mult_utils.h>

// function defined in ntt_scala.c
ring_t getRing(int degree);
void poly_fast_ntt_transform(ntt_t* dst, polynomial_t src, ring_t ring, int rootOfUnity);
void poly_fast_inv_ntt_tranform(polynomial_t* dst, ntt_t src, ring_t ring); 
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
        vec_acc = __riscv_vrem_vx_i32m8(vec_acc, dst->modulo, vl);
        // storing results
        __riscv_vse32_v_i32m8(dst_coeffs, vec_acc, vl);
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

    poly_fast_ntt_transform(&ntt_lhs, lhs, ring, ring.rootOfUnity);
    poly_fast_ntt_transform(&ntt_lhs_times_rhs, rhs, ring, ring.rootOfUnity);

    rvv_ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    poly_fast_inv_ntt_tranform(dst, ntt_lhs_times_rhs, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

/** Benchmark wrapper */
poly_mult_bench_result_t poly_mult_ntt_rvv_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv, golden);
}