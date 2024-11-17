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


/** Compute n-element NTT, assuming level @p level
 *
 *
 * @param[out] dst destination buffer for NTT transform result
 * @param[in] inputs coefficients (must contain @p n elements)
 * @param n number of coefficients
 * @param level NTT level (start from 0)
 * @param rootPowers 2D array of pre-computed root of unit powers rootPowers[level][i] = (rootOfUnit ^ level) ^ i
*/
void rvv_ntt_transform_helper(ntt_t* dst, int* coeffs, int n, int level, int rootPowers[8][64]) {
    const size_t coeffWidth = sizeof(coeffs[0]);

    // TODO: optimize for n > 1 (target is to optimize as soon as all coefficients fit in
    //       a vector register group)
    if (n == 1) {
        dst->coeffs[0] = coeffs[0];
        return;
    }

    size_t avl = n / 2; // half of n odd/even coefficients
    ntt_t ntt_even = allocate_poly(n / 2 - 1, dst->modulo);
    ntt_t ntt_odd = allocate_poly(n / 2 - 1, dst->modulo);
    int* even_coeffs = ntt_even.coeffs;
    int* odd_coeffs = ntt_odd.coeffs; // coeffs + n / 2;
    for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, coeffs += 2*vl)
    {
        vl = __riscv_vsetvl_e32m8(avl);
        // splitting even and odd coefficients using strided load
        vint32m8_t vec_even_coeffs = __riscv_vlse32_v_i32m8((int*) coeffs, 2 * sizeof(coeffs[0]), vl);
        vint32m8_t vec_odd_coeffs = __riscv_vlse32_v_i32m8((int*) (coeffs + 1), 2 * sizeof(coeffs[0]), vl);

        __riscv_vse32_v_i32m8((int*) even_coeffs, vec_even_coeffs, vl);
        __riscv_vse32_v_i32m8((int*) odd_coeffs, vec_odd_coeffs, vl);
    }

    // NTT recursion
    ntt_t dst_even = {.degree = (n / 2 - 1), .modulo = dst->modulo, coeffs = dst->coeffs, .coeffSize = (n/2)};
    ntt_t dst_odd = {.degree = (n / 2 - 1), .modulo = dst->modulo, coeffs = (dst->coeffs + n / 2), .coeffSize = (n/2)};
    rvv_ntt_transform_helper(&dst_even, ntt_even.coeffs, n / 2, level + 1, rootPowers);
    rvv_ntt_transform_helper(&dst_odd, ntt_odd.coeffs, n / 2, level + 1, rootPowers);

    free(ntt_even.coeffs);
    free(ntt_odd.coeffs);

    // final stage of the butterfly 
    avl = n / 2;
    assert(avl <= 64); // rootPowers[level] is at most a 64-element array
    even_coeffs = dst->coeffs;
    odd_coeffs = dst->coeffs + n / 2;
    int* twiddleFactor = rootPowers[level];
    for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, twiddleFactor += vl)
    {
        vl = __riscv_vsetvl_e32m8(avl);
        vint32m8_t vec_even_coeffs = __riscv_vle32_v_i32m8((int*) even_coeffs, vl);
        vint32m8_t vec_odd_coeffs = __riscv_vle32_v_i32m8((int*) odd_coeffs, vl);

        vint32m8_t vec_twiddleFactor = __riscv_vle32_v_i32m8((int*) twiddleFactor, vl);

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

#define LMUL 4
#define E32_MASK 8 // 32 / LMUL

#define BUILD_LMUL_ARGS(x) x, LMUL
#define BUILD_LMUL_REVARGS(x) LMUL, x
#define CONCAT(x, y) x ## y
#define CONCAT3(x, y, z) x ## y ## z
#define CONCAT_FWD(ARGS) CONCAT(ARGS)
#define CONCAT3_FWD(ARGS) CONCAT3(ARGS)
#define PREPEND_LMUL(x) CONCAT_FWD(BUILD_LMUL_ARGS(x))
#define APPEND_LMUL(x) CONCAT(BUILD_LMUL_REV_ARGS(x))
#define PREPEND_MASK_E32(x) CONCAT_FWD(BUILD_E32_MASK_ARGS(x))
#define BUILD_E32_MASK_ARGS(x) x, E32_MASK


#define _FUNC_LMUL(CMD, SUFFIX) CMD ## SUFFIX
#define _TYPE_LMUL(FMT, SUFFIX) FMT ## SUFFIX ## _t

#define _FUNC_LMUL_ARG1(ARGS) _FUNC_LMUL(ARGS)
#define _TYPE_LMUL_ARG1(ARGS) _TYPE_LMUL(ARGS)

#define FUNC_LMUL(CMD) _FUNC_LMUL_ARG1(BUILD_PARAMS(CMD))
#define FUNC_LMUL_MASKED(CMD) _FUNC_LMUL_ARG1(BUILD_PARAMS_LMUL_MASKED(CMD))
#define TYPE_LMUL(FMT) _TYPE_LMUL_ARG1(BUILD_PARAMS(FMT))
#define MASK_TYPE_E32(FMT) _TYPE_LMUL_ARG1(BUILD_PARAMS_MASK_E32_TYPE(FMT))
#define MASK_FUNC_E32(CMD) _FUNC_LMUL_ARG1(BUILD_PARAMS_MASK_E32_FUNC(CMD))

#define BUILD_PARAMS(FMT) FMT, PREPEND_LMUL(m)
#define BUILD_PARAMS_MASK_E32_TYPE(FMT) FMT, E32_MASK
#define BUILD_PARAMS_MASK_E32_FUNC(CMD) CMD, E32_MASK
#define BUILD_PARAMS_LMUL_MASKED(FMT) FMT, LMUL_MASKED_SUFFIX
#define LMUL_MASKED_SUFFIX CONCAT3_FWD(BUILD_LMUL_MASKED_SUFFIX_ARGS)
#define BUILD_LMUL_MASKED_SUFFIX_ARGS m, LMUL, _mu

// Select between strided loads and vcompress (+ unit-stride load) to perform the odd/even split
#ifndef USE_STRIDED_LOAD
#define USE_STRIDED_LOAD 0
#endif // defined(USE_STRIDED_LOAD)

// Select between implementation using a pre-computer array of root powers or not
#ifndef USE_PRECOMPUTED_ROOT_POWERS
#define USE_PRECOMPUTED_ROOT_POWERS 1
#endif // defined(USE_PRECOMPUTED_ROOT_POWERS)

/** Compute n-element NTT, assuming level @p level
 *
 *
 * @param[out] dst destination buffer for NTT transform result
 * @param[in] inputs coefficients (must contain @p n elements)
 * @param n number of coefficients
 * @param level NTT level (start from 0)
 * @param rootPowers 2D array of pre-computed root of unit powers rootPowers[level][i] = (rootOfUnit ^ level) ^ i
*/
void rvv_ntt_transform_fast_helper(ntt_t* dst, int* coeffs, int _n, int level, int rootPowers[8][64]) {
    const size_t coeffWidth = sizeof(coeffs[0]);

    size_t vlmax = FUNC_LMUL(__riscv_vsetvlmax_e32)();

    assert(_n > 1);

    // A/B buffering for coefficients
    int* coeffs_a = coeffs;
    int* coeffs_b = dst->coeffs;

    // iteration over the local number of coefficients
    int local_level;
    int n;
    // masks used for the odd/even split using vcompress instructions
    vint8m1_t mask_even_i8 = __riscv_vmv_v_x_i8m1(0x55, __riscv_vsetvlmax_e32m1());
    vint8m1_t mask_odd_i8 = __riscv_vmv_v_x_i8m1(0xAA, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_even_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_even_i8);
    MASK_TYPE_E32(vbool) mask_odd_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_odd_i8);

    for (local_level = level, n = _n; n > 4; n = n / 2, local_level++) {
        const int m = 1 << local_level;
        const int half_n = n / 2;

        for (int j = 0, coeffs_a_offset = 0; j < m; j++) {

            if (USE_STRIDED_LOAD) {
                size_t avl = half_n; // half of n odd/even coefficients
                int* even_coeffs = coeffs_b + 2 * j * half_n;
                int* odd_coeffs = even_coeffs + half_n;
                for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, coeffs_a_offset += 2*vl)
                {
                    vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
                    // splitting even and odd coefficients using strided load
                    TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vlse32_v_i32)((int*) coeffs_a + coeffs_a_offset, 2 * sizeof(coeffs[0]), vl);
                    TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vlse32_v_i32)((int*) (coeffs_a + coeffs_a_offset+ 1), 2 * sizeof(coeffs[0]), vl);

                    FUNC_LMUL(__riscv_vse32_v_i32)((int*) even_coeffs, vec_even_coeffs, vl);
                    FUNC_LMUL(__riscv_vse32_v_i32)((int*) odd_coeffs, vec_odd_coeffs, vl);
                } 
            } else {
                size_t avl = n; // half of n odd/even coefficients
                int* even_coeffs = coeffs_b + 2 * j * half_n;
                int* odd_coeffs = even_coeffs + half_n;
                for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl / 2, odd_coeffs += vl / 2, coeffs_a_offset += vl)
                {
                    // using a unit-stride load and a pair of vcompress for the even/odd split
                    vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
                    TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) coeffs_a + coeffs_a_offset, vl);
                    // 2^level coefficients
                    // mask type is bool<n> where n in EEW / EMUL (in our case 32 / 8 = 4)
                    TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vcompress_vm_i32)(vec_coeffs, mask_even_b4, vl);
                    TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vcompress_vm_i32)(vec_coeffs, mask_odd_b4, vl);

                    FUNC_LMUL(__riscv_vse32_v_i32)((int*) even_coeffs, vec_even_coeffs, vl / 2);
                    FUNC_LMUL(__riscv_vse32_v_i32)((int*) odd_coeffs, vec_odd_coeffs, vl / 2);
                }
            }
        }
        // swapping A/B buffers
        int* tmp = coeffs_a;
        coeffs_a = coeffs_b;
        coeffs_b = tmp;
    }

    // final levels
    assert(n == 4);
    vint8m1_t mask_id_i8 = __riscv_vmv_v_x_i8m1(0x99, __riscv_vsetvlmax_e32m1());
    vint8m1_t mask_up_i8 = __riscv_vmv_v_x_i8m1(0x44, __riscv_vsetvlmax_e32m1());
    vint8m1_t mask_down_i8 = __riscv_vmv_v_x_i8m1(0x22, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_id_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_i8);
    MASK_TYPE_E32(vbool) mask_up_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_i8);
    MASK_TYPE_E32(vbool) mask_down_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_down_i8);

    // 4-element butterfly unrolled across the full row dimension
    // assuming 4-elements fit into the vector register group, no vslideup nor down actually
    // move any active element outside of the vector register group
    for (; n > 2; n = n / 2, local_level++) {
        const int m = 1 << local_level;

        size_t avl = _n;
        int* coeffs_addr = coeffs_a;
        int* dst_coeffs = dst->coeffs; // making sure results are stored in the destination buffer
        for (size_t vl; avl > 0; avl -= vl, coeffs_addr += vl, dst_coeffs += vl)
        {
            // using a unit-stride load and a pair of vslide(up/down) to swap each pair of middle
            // coefficients within each group of 4 elements.
            vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
            TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) coeffs_addr, vl);
            TYPE_LMUL(vint32) res_coeffs = vec_coeffs; // a copy is needed to avoid overwriting the middle coefficients
            res_coeffs = FUNC_LMUL_MASKED(__riscv_vslideup_vx_i32)(mask_up_b4, res_coeffs, vec_coeffs, 1, vl);
            res_coeffs = FUNC_LMUL_MASKED(__riscv_vslidedown_vx_i32)(mask_down_b4, res_coeffs, vec_coeffs, 1, vl);

            FUNC_LMUL(__riscv_vse32_v_i32)((int*) dst_coeffs, res_coeffs, vl);
        }
    }

    for (n=2, local_level = local_level; local_level >= 0; n = 2 * n, local_level--) {
        const int m = 1 << local_level;
        const int half_n = n / 2;

        for (int j = 0; j < m; j++) {
            // final stage of the butterfly 
            size_t avl = half_n;
            assert(avl <= 64); // rootPowers[level] is at most a 64-element array
            int* even_coeffs = coeffs_a + 2 * j * half_n;
            int* odd_coeffs = even_coeffs + half_n;
            int* twiddleFactor = rootPowers[local_level];
            for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, twiddleFactor += vl)
            {
                vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
                TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) even_coeffs, vl);
                TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) odd_coeffs, vl);

                TYPE_LMUL(vint32) vec_twiddleFactor = FUNC_LMUL(__riscv_vle32_v_i32)((int*) twiddleFactor, vl);

                // TODO: consider using a vectorized version of Barret's reduction algorithm
                TYPE_LMUL(vint32) vec_odd_results = FUNC_LMUL(__riscv_vmul_vv_i32)(vec_odd_coeffs, vec_twiddleFactor, vl);
                TYPE_LMUL(vint32) vec_even_results = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
                // even results
                vec_even_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_even_results, dst->modulo, vl);
                FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);

                // odd results
                vec_odd_results = FUNC_LMUL(__riscv_vsub_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
                vec_odd_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_odd_results, dst->modulo, vl);
                FUNC_LMUL(__riscv_vse32_v_i32)(odd_coeffs, vec_odd_results, vl);
            }
        } 
    }
}

void initRootPowerTable(ring_t ring, int rootPowers[8][64], int rootOfUnity);

void rvv_ntt_transform(ntt_t* dst, int* coeffs, ring_t ring, int degree, int rootOfUnity) {
    const int n = degree + 1;
    assert(n == 128);
    int rootPowers[8][64];
    initRootPowerTable(ring, rootPowers, rootOfUnity);

    rvv_ntt_transform_helper(dst, coeffs, n, 0, rootPowers);
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


void poly_mult_ntt_rvv_v2(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(dst->degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    if (USE_PRECOMPUTED_ROOT_POWERS) {
        // using pre-computed root powers
        rvv_ntt_transform_helper(&ntt_lhs, lhs.coeffs, lhs.degree + 1, 0, ringPowers[0]);
        rvv_ntt_transform_helper(&ntt_lhs_times_rhs, rhs.coeffs, rhs.degree + 1, 0, ringPowers[0]);
    } else {
        rvv_ntt_transform(&ntt_lhs, lhs.coeffs, ring, lhs.degree, ring.rootOfUnity);
        rvv_ntt_transform(&ntt_lhs_times_rhs, rhs.coeffs, ring, rhs.degree, ring.rootOfUnity);
    }

    // element-size multiplication using RVV
    rvv_ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    if (1) {
        // using pre-computed inverse root powers
        rvv_ntt_transform_helper(dst, ntt_lhs_times_rhs.coeffs, dst->degree + 1, 0, ringInvPowers[0]);
    } else {
        rvv_ntt_transform(dst, ntt_lhs_times_rhs.coeffs, ring, ntt_lhs_times_rhs.degree, ring.invRootOfUnity);
    }
    // division by the degree
    rvv_ntt_degree_scaling(dst, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

/** Benchmark wrapper */
poly_mult_bench_result_t poly_mult_ntt_rvv_v2_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_v2, golden);
}


void poly_mult_ntt_rvv_v3(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(dst->degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    rvv_ntt_transform_fast_helper(&ntt_lhs, lhs.coeffs, lhs.degree + 1, 0, ringPowers[0]);
    rvv_ntt_transform_fast_helper(&ntt_lhs_times_rhs, rhs.coeffs, rhs.degree + 1, 0, ringPowers[0]);

    // element-size multiplication using RVV
    rvv_ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    rvv_ntt_transform_fast_helper(dst, ntt_lhs_times_rhs.coeffs, dst->degree + 1, 0, ringInvPowers[0]);

    // division by the degree
    rvv_ntt_degree_scaling(dst, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

/** Benchmark wrapper */
poly_mult_bench_result_t poly_mult_ntt_rvv_v3_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_v3, golden);
}