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

#ifndef LMUL
#define LMUL 2
#define WLMUL 4 // 2 * LMUL (for widening)
#define E32_MASK 16 // 32 / LMUL
#endif

#ifndef E32_MASK
#error "E32_MASK (32 / LMUL) must be defined"
#endif

#ifndef WLMUL
#error "WLMUL (2 * LMUL) must be defined"
#endif

#define BUILD_LMUL_ARGS(x) x, LMUL
#define BUILD_WLMUL_ARGS(x) x, WLMUL
#define CONCAT(x, y) x ## y
#define CONCAT3(x, y, z) x ## y ## z
#define CONCAT_FWD(ARGS) CONCAT(ARGS)
#define CONCAT3_FWD(ARGS) CONCAT3(ARGS)
#define PREPEND_LMUL(x) CONCAT_FWD(BUILD_LMUL_ARGS(x))
#define PREPEND_WLMUL(x) CONCAT_FWD(BUILD_WLMUL_ARGS(x))
#define APPEND_LMUL(x) CONCAT(BUILD_LMUL_REV_ARGS(x))


#define _FUNC_LMUL(CMD, SUFFIX) CMD ## SUFFIX
#define _TYPE_LMUL(FMT, SUFFIX) FMT ## SUFFIX ## _t
#define _FUNC_LMUL_2PART(CMD0, SUFFIX0, CMD1, SUFFIX1) CMD0 ## SUFFIX0 ## CMD1 ## SUFFIX1

#define _FUNC_LMUL_ARG1(ARGS) _FUNC_LMUL(ARGS)
#define _TYPE_LMUL_ARG1(ARGS) _TYPE_LMUL(ARGS)
#define _FUNC_LMUL_2PART_ARG1(ARGS) _FUNC_LMUL_2PART(ARGS)

#define FUNC_LMUL(CMD) _FUNC_LMUL_ARG1(BUILD_PARAMS(CMD))
#define WFUNC_LMUL(CMD) _FUNC_LMUL_ARG1(BUILD_WPARAMS(CMD))
#define FUNC_LMUL_2PART(CMD0,CMD1) _FUNC_LMUL_2PART_ARG1(BUILD_PARAMS_2PART(CMD0,CMD1))
#define FUNC_LMUL_MASKED(CMD) _FUNC_LMUL_ARG1(BUILD_PARAMS_LMUL_MASKED(CMD))
#define TYPE_LMUL(FMT) _TYPE_LMUL_ARG1(BUILD_PARAMS(FMT))
#define WTYPE_LMUL(FMT) _TYPE_LMUL_ARG1(BUILD_WPARAMS(FMT))
#define MASK_TYPE_E32(FMT) _TYPE_LMUL_ARG1(BUILD_PARAMS_MASK_E32(FMT))
#define MASK_FUNC_E32(CMD) _FUNC_LMUL_ARG1(BUILD_PARAMS_MASK_E32(CMD))
#define MASK_LMUL_FUNC_E32(CMD) _FUNC_LMUL_2PART_ARG1(BUILD_PARAMS_MASK_E32_LMUL_FUNC(CMD))

#define BUILD_PARAMS(FMT) FMT, PREPEND_LMUL(m)
#define BUILD_WPARAMS(FMT) FMT, PREPEND_WLMUL(m)
#define BUILD_PARAMS_2PART(CMD0,CMD1) CMD0, PREPEND_LMUL(m), CMD1, PREPEND_LMUL(m)
#define BUILD_PARAMS_MASK_E32(FMT) FMT, E32_MASK
#define BUILD_PARAMS_MASK_E32_LMUL_FUNC(CMD) CMD, PREPEND_LMUL(m), _b, E32_MASK
#define BUILD_PARAMS_LMUL_MASKED(FMT) FMT, LMUL_MASKED_SUFFIX
#define LMUL_MASKED_SUFFIX CONCAT3_FWD(BUILD_LMUL_MASKED_SUFFIX_ARGS)
#define BUILD_LMUL_MASKED_SUFFIX_ARGS m, LMUL, _mu

// Select between implementation using a pre-computer array of root powers or not
#ifndef USE_PRECOMPUTED_ROOT_POWERS
#define USE_PRECOMPUTED_ROOT_POWERS 1
#endif // defined(USE_PRECOMPUTED_ROOT_POWERS)

// Coefficient indices generated with script poly_mult_emulation.py
// NTT forward pass indices for 8 coeffs:
const int32_t ntt_coeff_indices_8[] = {
   0, 16, 8, 24, 4, 20, 12, 28 
};
// NTT forward pass indices for 16 coeffs:
const int32_t ntt_coeff_indices_16[] = {
   0, 32, 16, 48, 8, 40, 24, 56,
  4, 36, 20, 52, 12, 44, 28, 60 
};
// NTT forward pass indices for 32 coeffs:
const int32_t ntt_coeff_indices_32[] = {
   0, 64, 32, 96, 16, 80, 48, 112,
  8, 72, 40, 104, 24, 88, 56, 120,
  4, 68, 36, 100, 20, 84, 52, 116,
  12, 76, 44, 108, 28, 92, 60, 124 
};
// NTT forward pass indices for 64 coeffs:
const int32_t ntt_coeff_indices_64[] = {
   0, 128, 64, 192, 32, 160, 96, 224,
  16, 144, 80, 208, 48, 176, 112, 240,
  8, 136, 72, 200, 40, 168, 104, 232,
  24, 152, 88, 216, 56, 184, 120, 248,
  4, 132, 68, 196, 36, 164, 100, 228,
  20, 148, 84, 212, 52, 180, 116, 244,
  12, 140, 76, 204, 44, 172, 108, 236,
  28, 156, 92, 220, 60, 188, 124, 252 
};
// NTT forward pass indices for 128 coeffs:
const uint32_t ntt_coeff_indices_128[] = {
   0, 256, 128, 384, 64, 320, 192, 448,
  32, 288, 160, 416, 96, 352, 224, 480,
  16, 272, 144, 400, 80, 336, 208, 464,
  48, 304, 176, 432, 112, 368, 240, 496,
  8, 264, 136, 392, 72, 328, 200, 456,
  40, 296, 168, 424, 104, 360, 232, 488,
  24, 280, 152, 408, 88, 344, 216, 472,
  56, 312, 184, 440, 120, 376, 248, 504,
  4, 260, 132, 388, 68, 324, 196, 452,
  36, 292, 164, 420, 100, 356, 228, 484,
  20, 276, 148, 404, 84, 340, 212, 468,
  52, 308, 180, 436, 116, 372, 244, 500,
  12, 268, 140, 396, 76, 332, 204, 460,
  44, 300, 172, 428, 108, 364, 236, 492,
  28, 284, 156, 412, 92, 348, 220, 476,
  60, 316, 188, 444, 124, 380, 252, 508 
};


// Function to display each element of a vint32m1_t vector
void display_vint32(TYPE_LMUL(vint32) vec, size_t vl) {
    int32_t elements[vl];
    FUNC_LMUL(__riscv_vse32_v_i32)(elements, vec, vl);

    printf("Vector elements: ");
    for (size_t i = 0; i < vl; i++) {
        printf("%d ", elements[i]);
    }
    printf("\n");
}

/** RVV-based implementation of Barrett's method for modulo reduction (modulo 3329)
 *
 * @param v input vector
 * @param vl vector length (number of element modulo reduced)
 * @return remainder
 */
static inline TYPE_LMUL(vint32) rvv_barrett_reduction(TYPE_LMUL(vint32) v, size_t vl) {
    // int mul = (lhs.coeffs[d] * rhs.coeffs[d]);
    // int tmp = mul - ((((int64_t) mul * 5039LL) >> 24)) * 3329;
    // dst->coeffs[d] = tmp >= 3329 ? tmp - 3329 : tmp; 
    // v is the result of a multiplication between two coefficients (each up to 13-bit wide)
    // so 64-bit dynamic is required when multiplying it with 5039 (13-bit unsigned)
    WTYPE_LMUL(vint64) vec_wide_results = WFUNC_LMUL(__riscv_vwmul_vx_i64)(v, 5039, vl);
    TYPE_LMUL(vint32) vec_tmp_results = FUNC_LMUL(__riscv_vnsra_wx_i32)(vec_wide_results, 24, vl);
    v = FUNC_LMUL(__riscv_vnmsac_vx_i32)(v, 3329, vec_tmp_results, vl);
    MASK_TYPE_E32(vbool) cmp_mask = MASK_LMUL_FUNC_E32(__riscv_vmsge_vx_i32)(v, 3329, vl);
    v = FUNC_LMUL_MASKED(__riscv_vadd_vx_i32)(cmp_mask, v, v, -3329, vl);

    return v;
}

static inline TYPE_LMUL(vint32) rvv_ntt_butterfly(TYPE_LMUL(vint32) vec_coeffs, int n, int modulo, TYPE_LMUL(vint32) vec_twiddleFactor, MASK_TYPE_E32(vbool) mask_up_b4, size_t vl) {
    vec_coeffs = FUNC_LMUL_MASKED(__riscv_vmul_vv_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vec_twiddleFactor, vl);
    
    // swapping odd/even pairs of coefficients
    TYPE_LMUL(vint32) vec_swapped_coeffs = FUNC_LMUL(__riscv_vslidedown_vx_i32)(vec_coeffs, n / 2, vl);
    vec_swapped_coeffs = FUNC_LMUL_MASKED(__riscv_vslideup_vx_i32)(mask_up_b4, vec_swapped_coeffs, vec_coeffs, n / 2, vl);

    vec_coeffs = FUNC_LMUL_MASKED(__riscv_vneg_v_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vl);
    vec_coeffs = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_coeffs, vec_swapped_coeffs, vl);

    // TODO: optimize modulo reduction
#if 0
    vec_coeffs = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_coeffs, modulo, vl);
#else
    vec_coeffs = rvv_barrett_reduction(vec_coeffs, vl);
#endif

    return vec_coeffs;
}

typedef struct {
    int _USE_STRIDED_LOAD;
    int _USE_INDEXED_LOAD;
    int _FINAL_N;
    int _BARRETT_RED;
    int _UNROLL_STOP;
    int _FUSED_BUTTERFLY;
} ntt_params_t;

/** Compute n-element NTT, assuming level @p level
 *  WARNING: non-indexed variant is destructive as coeffs is used as a temporary buffer and overwritten
 *
 * @param[out] dst destination buffer for NTT transform result
 * @param[in] inputs coefficients (must contain @p n elements)
 * @param n number of coefficients
 * @param level NTT level (start from 0)
 * @param rootPowers 2D array of pre-computed root of unit powers rootPowers[level][i] = (rootOfUnit ^ level) ^ i
*/
void rvv_ntt_transform_fast_helper(ntt_t* dst, int* coeffs, int _n, int level, int rootPowers[8][64], ntt_params_t params) {
    const size_t coeffWidth = sizeof(coeffs[0]);

    assert(_n > 1);

    // A/B buffering for coefficients
    // TODO/FIXM: non-indexed variant overwrite the coeffs input array
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

    if (params._USE_INDEXED_LOAD)
    {
        size_t avl = _n; // half of n odd/even coefficients
        const unsigned int* coeffs_index = ntt_coeff_indices_128;
        int* dst_coeffs = dst->coeffs;
        for (size_t vl; avl > 0; avl -= vl, coeffs_index += vl, dst_coeffs += vl)
        {
            vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);

            TYPE_LMUL(vuint32) vec_indices = FUNC_LMUL(__riscv_vle32_v_u32)(coeffs_index, vl);
            TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vluxei32_v_i32)((const int*) coeffs_a, vec_indices, vl);

            FUNC_LMUL(__riscv_vse32_v_i32)((int*) dst_coeffs, vec_coeffs, vl);
        }

        local_level = 6;
        n = 2;

    } else {
        // initial stage(s) of the explicit coefficient permutations
        for (local_level = level, n = _n; n > params._FINAL_N; n = n / 2, local_level++) {
            const int m = 1 << local_level;
            const int half_n = n / 2;

            for (int j = 0, coeffs_a_offset = 0; j < m; j++) {

                if (params._USE_STRIDED_LOAD) {
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
        if (n == 4) { // allowing the last level to be split out of the previous 2D loop and optimized independently
            assert(n == 4);
            vint8m1_t mask_up_i8 = __riscv_vmv_v_x_i8m1(0x44, __riscv_vsetvlmax_e32m1());
            vint8m1_t mask_down_i8 = __riscv_vmv_v_x_i8m1(0x22, __riscv_vsetvlmax_e32m1());
            MASK_TYPE_E32(vbool) mask_up_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_i8);
            MASK_TYPE_E32(vbool) mask_down_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_down_i8);

            // final stage(s) of the explicit coefficient permutations
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
        } else {
            if (coeffs_a != dst->coeffs) {
                // copying the result back to the destination
                memcpy(dst->coeffs, coeffs_a, _n * sizeof(int));
                // TODO: this can certainly be optimized by having different source/destination pointers
                //       for the first local level of the next loop if the buffer differs
            }
        }
    }
    // reconstruction stage input should be equal to the destination buffer
    // (a copy was inserted to ensure this is true)
    coeffs_a = dst->coeffs;

    // condition under which butterfly level 4 can be optimized (required a wide enought vector register group)
    const int cond_lvl4 = (params._UNROLL_STOP < 4 && 8 <= FUNC_LMUL(__riscv_vsetvlmax_e32)());
    // mask for 4th level of butterfly (used for vslideup to swap 4-elt sub-group in 8-elt group)
    vint8m1_t mask_up_lvl4_i8 = __riscv_vmv_v_x_i8m1(0xf0, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_up_lvl4_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_lvl4_i8);
    // expect rootPowers built with replicate parameter set so the 2-element pattern is repeated across the vector register group
    // local_level is 4 (used to index rootPowers)
    TYPE_LMUL(vint32) vec_twiddleFactor_lvl4 = FUNC_LMUL(__riscv_vle32_v_i32)((int32_t*) rootPowers[4], FUNC_LMUL(__riscv_vsetvlmax_e32)());

    const int cond_lvl5 = (params._UNROLL_STOP < 5 && 4 <= FUNC_LMUL(__riscv_vsetvlmax_e32)());
    vint8m1_t mask_up_lvl5_i8 = __riscv_vmv_v_x_i8m1(0xcc, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_up_lvl5_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_lvl5_i8);

    // expect rootPowers built with replicate parameter set so the 2-element pattern is repeated across the vector register group
    TYPE_LMUL(vint32) vec_twiddleFactor_lvl5 = FUNC_LMUL(__riscv_vle32_v_i32)((int32_t*) rootPowers[5], FUNC_LMUL(__riscv_vsetvlmax_e32)());

    assert(n == 2);
    // optimize n=2 case (vslide(up/down) which 0x55 and 0xAA masks)
    if (params._UNROLL_STOP < 6) {
        assert (n == 2 && local_level == 6);
        size_t avl = _n;
        int* coeffs_addr = coeffs_a;

        vint8m1_t mask_up_i8 = __riscv_vmv_v_x_i8m1(0xaa, __riscv_vsetvlmax_e32m1());
        MASK_TYPE_E32(vbool) mask_up_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_i8);

        // first level, no need to multiply by any twiddle factor
        // because this level the factor is twiddleFactor^0 = 1
        for (size_t vl; avl > 0; avl -= vl, coeffs_addr += vl)
        {
            vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);

            TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) coeffs_addr, vl);
            
            // swapping odd/even pairs of coefficients
            TYPE_LMUL(vint32) vec_swapped_coeffs = FUNC_LMUL(__riscv_vslidedown_vx_i32)(vec_coeffs, 1, vl);
            vec_swapped_coeffs = FUNC_LMUL_MASKED(__riscv_vslideup_vx_i32)(mask_up_b4, vec_swapped_coeffs, vec_coeffs, 1, vl);

            vec_coeffs = FUNC_LMUL_MASKED(__riscv_vneg_v_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vl);
            vec_coeffs = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_coeffs, vec_swapped_coeffs, vl);

            if (params._FUSED_BUTTERFLY && cond_lvl5) {
                // fusing optimization for the next level (4-elt group case)
                vec_coeffs = rvv_ntt_butterfly(vec_coeffs, 2 * n, dst->modulo, vec_twiddleFactor_lvl5, mask_up_lvl5_b4, vl);
            }

            if (params._FUSED_BUTTERFLY && cond_lvl4) {
                // fusing optimization for the next level (8-elt group case)
                vec_coeffs = rvv_ntt_butterfly(vec_coeffs, 4 * n, dst->modulo, vec_twiddleFactor_lvl4, mask_up_lvl4_b4, vl);
            }

            // Note: no modulo reduction is performed
            FUNC_LMUL(__riscv_vse32_v_i32)(coeffs_addr, vec_coeffs, vl);
        }
        n = 2 * n;
        local_level--;

        if (params._FUSED_BUTTERFLY && cond_lvl5) {
            assert(n == 4 && local_level == 5);
            assert(n <= FUNC_LMUL(__riscv_vsetvlmax_e32)()); // the group number of elements n must fit in vlmax for current LMUL
            // extra level fused
            n = 2 * n;
            local_level--;
        }

        if (params._FUSED_BUTTERFLY && cond_lvl4) {
            assert(n == 8 && local_level == 4);
            assert(n <= FUNC_LMUL(__riscv_vsetvlmax_e32)()); // the group number of elements n must fit in vlmax for current LMUL
            // extra level fused
            n = 2 * n;
            local_level--;
        }
    }

    // optimize n=4 case (vslide(up/down) which 0x33 and 0xcc masks)
    if (!params._FUSED_BUTTERFLY && cond_lvl5) {
        assert (n == 4 && local_level == 5);
        size_t avl = _n;
        int* coeffs_addr = coeffs_a;

        for (size_t vl; avl > 0; avl -= vl, coeffs_addr += vl)
        {
            vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);

            TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) coeffs_addr, vl);

            vec_coeffs = rvv_ntt_butterfly(vec_coeffs, n, dst->modulo, vec_twiddleFactor_lvl5, mask_up_lvl5_b4, vl);

            if (params._FUSED_BUTTERFLY && cond_lvl4) {
                // fusing optimization for the next level (8-elt group case)
                vec_coeffs = rvv_ntt_butterfly(vec_coeffs, 2 * n, dst->modulo, vec_twiddleFactor_lvl4, mask_up_lvl4_b4, vl);
            }

            FUNC_LMUL(__riscv_vse32_v_i32)(coeffs_addr, vec_coeffs, vl);
        }
        n = 2 * n;
        local_level--;

        if (params._FUSED_BUTTERFLY && cond_lvl4) {
            assert(n == 8 && local_level == 4);
            assert(n <= FUNC_LMUL(__riscv_vsetvlmax_e32)()); // the group number of elements n must fit in vlmax for current LMUL
            // extra level fused
            n = 2 * n;
            local_level--;
        }
    }


    // optimize n=8 case (vslide(up/down) which 0x0f and 0xf0 masks)
    if (!params._FUSED_BUTTERFLY && cond_lvl4) {
        assert (n == 8 && local_level == 4);
        assert(n <= FUNC_LMUL(__riscv_vsetvlmax_e32)()); // the group number of elements n must fit in vlmax for current LMUL
        size_t avl = _n;
        int* coeffs_addr = coeffs_a;

        avl = _n;

        for (size_t vl; avl > 0; avl -= vl, coeffs_addr += vl)
        {
            vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
            TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) coeffs_addr, vl);

            vec_coeffs = rvv_ntt_butterfly(vec_coeffs, n, dst->modulo, vec_twiddleFactor_lvl4, mask_up_lvl4_b4, vl);
            FUNC_LMUL(__riscv_vse32_v_i32)(coeffs_addr, vec_coeffs, vl);
        }
        n = 2 * n;
        local_level--;
    }

    // disabling assert to allow skipping some optimization levels due to LMUL too small to fit
    // a local element group to perform swapping 
    // assert(local_level == params._UNROLL_STOP && n == 2 << (6 - params._UNROLL_STOP));
    for (; local_level >= 0; n = 2 * n, local_level--) {
        const int m = 1 << local_level;
        const int half_n = n / 2;

        for (int j = 0; j < m; j++) {
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

                TYPE_LMUL(vint32) vec_odd_results;
                TYPE_LMUL(vint32) vec_even_results;

                vec_odd_results = FUNC_LMUL(__riscv_vmul_vv_i32)(vec_odd_coeffs, vec_twiddleFactor, vl);
                vec_even_results = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
                vec_odd_results = FUNC_LMUL(__riscv_vsub_vv_i32)(vec_even_coeffs, vec_odd_results, vl);

                if (params._BARRETT_RED) {
                    vec_odd_results = rvv_barrett_reduction(vec_odd_results, vl);
                    vec_even_results = rvv_barrett_reduction(vec_even_results, vl);
                } else {
                    // even results
                    vec_even_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_even_results, dst->modulo, vl);

                    // odd results
                    vec_odd_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_odd_results, dst->modulo, vl);
                }

                FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);
                FUNC_LMUL(__riscv_vse32_v_i32)(odd_coeffs, vec_odd_results, vl);
            }
        } 
    }
}


void initRootPowerTable(ring_t ring, int rootPowers[8][64], int rootOfUnity);

/** Compute n-element NTT, assuming level @p level
 *  WARNING: non-indexed variant is destructive as coeffs is used as a temporary buffer and overwritten
 *
 * @param[out] dst destination buffer for NTT transform result
 * @param[in] inputs coefficients (must contain @p n elements)
 * @param n number of coefficients
 * @param level NTT level (start from 0)
 * @param rootPowers 2D array of pre-computed root of unit powers rootPowers[level][i] = (rootOfUnit ^ level) ^ i
*/
void rvv_ntt_transform_fastest_helper(ntt_t* dst, int* coeffs, int _n, int level, int rootPowers[8][64]) {
    const size_t coeffWidth = sizeof(coeffs[0]);

    assert(_n > 1);

    // iteration over the local number of coefficients
    int local_level = 6;
    int n = _n;
    // masks used for the odd/even split using vcompress instructions
    vint8m1_t mask_odd_i8 = __riscv_vmv_v_x_i8m1(0xAA, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_odd_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_odd_i8);

    {
        size_t avl = _n;
        const unsigned int* coeffs_index = ntt_coeff_indices_128;
        int* dst_coeffs = dst->coeffs;
        for (size_t vl; avl > 0; avl -= vl, coeffs_index += vl, dst_coeffs += vl)
        {
            vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);

            TYPE_LMUL(vuint32) vec_indices = FUNC_LMUL(__riscv_vle32_v_u32)(coeffs_index, vl);
            TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vluxei32_v_i32)((const int*) coeffs, vec_indices, vl);

            FUNC_LMUL(__riscv_vse32_v_i32)((int*) dst_coeffs, vec_coeffs, vl);
        }

        local_level = 6;
        n = 2;
    }

    return;
    // reconstruction stage input should be equal to the destination buffer
    // (a copy was inserted to ensure this is true)
    int* coeffs_a = dst->coeffs;

    vint8m1_t mask_up_lvl6_i8 = __riscv_vmv_v_x_i8m1(0xaa, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_up_lvl6_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_lvl6_i8);

    // condition under which butterfly level 4 can be optimized (required a wide enought vector register group)
    const int cond_lvl4 = (8 <= FUNC_LMUL(__riscv_vsetvlmax_e32)());
    // mask for 4th level of butterfly (used for vslideup to swap 4-elt sub-group in 8-elt group)
    vint8m1_t mask_up_lvl4_i8 = __riscv_vmv_v_x_i8m1(0xf0, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_up_lvl4_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_lvl4_i8);
    // expect rootPowers built with replicate parameter set so the 2-element pattern is repeated across the vector register group
    // local_level is 4 (used to index rootPowers)
    TYPE_LMUL(vint32) vec_twiddleFactor_lvl4 = FUNC_LMUL(__riscv_vle32_v_i32)((int32_t*) rootPowers[4], FUNC_LMUL(__riscv_vsetvlmax_e32)());

    const int cond_lvl5 = (4 <= FUNC_LMUL(__riscv_vsetvlmax_e32)());
    vint8m1_t mask_up_lvl5_i8 = __riscv_vmv_v_x_i8m1(0xcc, __riscv_vsetvlmax_e32m1());
    MASK_TYPE_E32(vbool) mask_up_lvl5_b4 = MASK_FUNC_E32(__riscv_vreinterpret_v_i8m1_b)(mask_up_lvl5_i8);

    // expect rootPowers built with replicate parameter set so the 2-element pattern is repeated across the vector register group
    TYPE_LMUL(vint32) vec_twiddleFactor_lvl5 = FUNC_LMUL(__riscv_vle32_v_i32)((int32_t*) rootPowers[5], FUNC_LMUL(__riscv_vsetvlmax_e32)());

    // optimize n=2, 4 and 8 cases  using masked vslide(up/down)
    {
        assert (n == 2 && local_level == 6);
        size_t avl = _n;
        int* coeffs_addr = coeffs_a;

        // first level, no need to multiply by any twiddle factor
        // because this level the factor is twiddleFactor^0 = 1
        for (size_t vl; avl > 0; avl -= vl, coeffs_addr += vl)
        {
            vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);

            TYPE_LMUL(vint32) vec_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) coeffs_addr, vl);
            
            // swapping odd/even pairs of coefficients
            TYPE_LMUL(vint32) vec_swapped_coeffs = FUNC_LMUL(__riscv_vslidedown_vx_i32)(vec_coeffs, 1, vl);
            vec_swapped_coeffs = FUNC_LMUL_MASKED(__riscv_vslideup_vx_i32)(mask_up_lvl6_b4, vec_swapped_coeffs, vec_coeffs, 1, vl);

            vec_coeffs = FUNC_LMUL_MASKED(__riscv_vneg_v_i32)(mask_up_lvl6_b4, vec_coeffs, vec_coeffs, vl);
            vec_coeffs = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_coeffs, vec_swapped_coeffs, vl);
            // Note: no modulo reduction is performed for level 6

            // fusing optimization for the next level (4-elt group case)
            vec_coeffs = rvv_ntt_butterfly(vec_coeffs, 4, dst->modulo, vec_twiddleFactor_lvl5, mask_up_lvl5_b4, vl);

            // fusing optimization for the next level (8-elt group case)
            vec_coeffs = rvv_ntt_butterfly(vec_coeffs, 8, dst->modulo, vec_twiddleFactor_lvl4, mask_up_lvl4_b4, vl);

            FUNC_LMUL(__riscv_vse32_v_i32)(coeffs_addr, vec_coeffs, vl);
        }
        n = 16;
        local_level = 3;

        assert(4 <= FUNC_LMUL(__riscv_vsetvlmax_e32)()); // the group number of elements n=4 for level 5 must fit in vlmax for current LMUL
        assert(8 <= FUNC_LMUL(__riscv_vsetvlmax_e32)()); // the group number of elements n=8 for level 4 must fit in vlmax for current LMUL
    }

    // disabling assert to allow skipping some optimization levels due to LMUL too small to fit
    // a local element group to perform swapping 
    for (; local_level >= 0; n = 2 * n, local_level--) {
        const int m = 1 << local_level;
        const int half_n = n / 2;

        for (int j = 0; j < m; j++) {
            size_t avl = half_n;
            assert(avl <= 64); // rootPowers[level] is at most a 64-element array
            int* even_coeffs = coeffs_a + 2 * j * half_n;
            int* odd_coeffs = even_coeffs + half_n;
            int* twiddleFactor = rootPowers[local_level];
            for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, twiddleFactor += vl)
            {
                vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
                // TODO even coefficients should be loaded after odd+twiddleFactor as they are needed last
                TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) odd_coeffs, vl);

                TYPE_LMUL(vint32) vec_twiddleFactor = FUNC_LMUL(__riscv_vle32_v_i32)((int*) twiddleFactor, vl);

                TYPE_LMUL(vint32) vec_odd_results = FUNC_LMUL(__riscv_vmul_vv_i32)(vec_odd_coeffs, vec_twiddleFactor, vl);
                TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) even_coeffs, vl);
                TYPE_LMUL(vint32) vec_even_results = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
                vec_odd_results = FUNC_LMUL(__riscv_vsub_vv_i32)(vec_even_coeffs, vec_odd_results, vl);

                if (1) {
                    vec_odd_results = rvv_barrett_reduction(vec_odd_results, vl);
                    vec_even_results = rvv_barrett_reduction(vec_even_results, vl);
                } else {
                    // even results
                    vec_even_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_even_results, dst->modulo, vl);

                    // odd results
                    vec_odd_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_odd_results, dst->modulo, vl);
                }

                FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);
                FUNC_LMUL(__riscv_vse32_v_i32)(odd_coeffs, vec_odd_results, vl);
            }
        } 
    }
}


void rvv_ntt_transform(ntt_t* dst, int* coeffs, ring_t ring, int degree, int rootOfUnity) {
    const int n = degree + 1;
    assert(n == 128);
    int rootPowers[8][64];
    initRootPowerTable(ring, rootPowers, rootOfUnity);

    rvv_ntt_transform_helper(dst, coeffs, n, 0, rootPowers);
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

    size_t avl = lhs.degree + 1;
    int* lhs_coeffs = lhs.coeffs;
    int* rhs_coeffs = rhs.coeffs;
    int* dst_coeffs = dst->coeffs;
    for (size_t vl; avl > 0; avl -= vl, lhs_coeffs += vl, rhs_coeffs += vl, dst_coeffs += vl)
    {
        // compute loop body vector length from avl (application vector length)
        vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
        // loading operands
        TYPE_LMUL(vint32) vec_src_lhs = FUNC_LMUL(__riscv_vle32_v_i32)(lhs_coeffs, vl);
        TYPE_LMUL(vint32) vec_src_rhs = FUNC_LMUL(__riscv_vle32_v_i32)(rhs_coeffs, vl);
        // modulo multiplication (eventually we will want to consider other techniques
        // than a naive remainder; e.g. Barret's reduction algorithm using a pre-computed
        // factor from the static modulo).
        TYPE_LMUL(vint32) vec_acc = FUNC_LMUL(__riscv_vmul_vv_i32)(vec_src_lhs, vec_src_rhs, vl);
#       if 0   
        vec_acc = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_acc, dst->modulo, vl);
#      else
        vec_acc = rvv_barrett_reduction(vec_acc, vl);
#      endif
        // storing results
        FUNC_LMUL(__riscv_vse32_v_i32)(dst_coeffs, vec_acc, vl);
    }
}

/** dividing each coefficient by the degree */
void rvv_ntt_degree_scaling(ntt_t* dst, ring_t ring) {
    size_t avl = dst->degree + 1;
    int* dst_coeffs = dst->coeffs;
    for (size_t vl; avl > 0; avl -= vl, dst_coeffs += vl)
    {
        // compute loop body vector length from avl (application vector length)
        vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
        // loading operands
        TYPE_LMUL(vint32) vec = FUNC_LMUL(__riscv_vle32_v_i32)(dst_coeffs, vl);
        // modulo multiplication (eventually we will want to consider other techniques
        // than a naive remainder; e.g. Barret's reduction algorithm using a pre-computed
        // factor from the static modulo).
        vec = FUNC_LMUL(__riscv_vmul_vx_i32)(vec, ring.invDegree, vl);
        // modulo reduction
#       if 0   
        vec = FUNC_LMUL(__riscv_vrem_vx_i32)(vec, ring.modulo, vl);
#      else
        vec = rvv_barrett_reduction(vec, vl);
#      endif
        // storing results
        FUNC_LMUL(__riscv_vse32_v_i32)(dst_coeffs, vec, vl);
    }
}

void poly_mult_ntt_rvv(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(lhs.degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    poly_fast_ntt_transform_helper(&ntt_lhs, lhs.coeffs, ring, lhs.degree, 1, 0, ringPowers[0]);
    poly_fast_ntt_transform_helper(&ntt_lhs_times_rhs, rhs.coeffs, ring, rhs.degree, 1, 0, ringPowers[0]);

    // element-size multiplication using RVV
    rvv_ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    poly_fast_ntt_transform_helper(dst, ntt_lhs_times_rhs.coeffs, ring, ntt_lhs_times_rhs.degree, 1, 0, ringInvPowers[0]);
    // division by the degree
    dst->degree = lhs.degree;
    rvv_ntt_degree_scaling(dst, ring);


    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

/** Benchmark wrapper */
poly_mult_bench_result_t poly_mult_ntt_rvv_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv, golden);
}


void poly_mult_ntt_rvv_recursive(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(lhs.degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(lhs.degree, 3329); 

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
        rvv_ntt_transform_helper(dst, ntt_lhs_times_rhs.coeffs, lhs.degree + 1, 0, ringInvPowers[0]);
    } else {
        rvv_ntt_transform(dst, ntt_lhs_times_rhs.coeffs, ring, ntt_lhs_times_rhs.degree, ring.invRootOfUnity);
    }
    // division by the degree
    dst->degree = lhs.degree;
    rvv_ntt_degree_scaling(dst, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

/** Benchmark wrapper */
poly_mult_bench_result_t poly_mult_ntt_rvv_recursive_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_recursive, golden);
}


void poly_mult_ntt_rvv_v3(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo, ntt_params_t params) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(lhs.degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(lhs.degree, 3329); 

    rvv_ntt_transform_fast_helper(&ntt_lhs, lhs.coeffs, lhs.degree + 1, 0, ringPowers[0], params);
    rvv_ntt_transform_fast_helper(&ntt_lhs_times_rhs, rhs.coeffs, rhs.degree + 1, 0, ringPowers[0], params);

    // element-size multiplication using RVV
    rvv_ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    rvv_ntt_transform_fast_helper(dst, ntt_lhs_times_rhs.coeffs, lhs.degree + 1, 0, ringInvPowers[0], params);

    // division by the degree
    dst->degree = lhs.degree;
    rvv_ntt_degree_scaling(dst, ring);

    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}


void poly_mult_ntt_rvv_fastest(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(lhs.degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(lhs.degree, 3329); 

    rvv_ntt_transform_fastest_helper(&ntt_lhs, lhs.coeffs, lhs.degree + 1, 0, ringPowers[0]);
    rvv_ntt_transform_fastest_helper(&ntt_lhs_times_rhs, rhs.coeffs, rhs.degree + 1, 0, ringPowers[0]);

    // element-size multiplication using RVV
    rvv_ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    rvv_ntt_transform_fastest_helper(dst, ntt_lhs_times_rhs.coeffs, lhs.degree + 1, 0, ringInvPowers[0]);

    // division by the degree
    dst->degree = lhs.degree;
    rvv_ntt_degree_scaling(dst, ring);

    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

void poly_mult_ntt_rvv_strided(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    ntt_params_t params_strided = {._USE_STRIDED_LOAD = 1, ._USE_INDEXED_LOAD = 0, ._FINAL_N = 4, ._BARRETT_RED = 0, ._UNROLL_STOP = 6, ._FUSED_BUTTERFLY = 0};
    poly_mult_ntt_rvv_v3(dst, lhs, rhs, modulo, params_strided);
}

void poly_mult_ntt_rvv_strided_barrett(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    ntt_params_t params_strided = {._USE_STRIDED_LOAD = 1, ._USE_INDEXED_LOAD = 0, ._FINAL_N = 4, ._BARRETT_RED = 1, ._UNROLL_STOP = 6, ._FUSED_BUTTERFLY = 0};
    poly_mult_ntt_rvv_v3(dst, lhs, rhs, modulo, params_strided);
}

void poly_mult_ntt_rvv_compressed(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    ntt_params_t params_compressed = {._USE_STRIDED_LOAD = 0, ._USE_INDEXED_LOAD = 0, ._FINAL_N = 4, ._BARRETT_RED = 0, ._UNROLL_STOP = 6, ._FUSED_BUTTERFLY = 0};
    poly_mult_ntt_rvv_v3(dst, lhs, rhs, modulo, params_compressed);
}

void poly_mult_ntt_rvv_indexed(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    ntt_params_t params_indexed = {._USE_STRIDED_LOAD = 0, ._USE_INDEXED_LOAD = 1, ._FINAL_N = 4, ._BARRETT_RED = 0, ._UNROLL_STOP = 6, ._FUSED_BUTTERFLY = 0};
    poly_mult_ntt_rvv_v3(dst, lhs, rhs, modulo, params_indexed);
}

void poly_mult_ntt_rvv_compressed_barrett(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    ntt_params_t params_compressed = {._USE_STRIDED_LOAD = 0, ._USE_INDEXED_LOAD = 0, ._FINAL_N = 4, ._BARRETT_RED = 1, ._UNROLL_STOP = 3, ._FUSED_BUTTERFLY = 1};
    poly_mult_ntt_rvv_v3(dst, lhs, rhs, modulo, params_compressed);
}

void poly_mult_ntt_rvv_indexed_barrett(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    ntt_params_t params_compressed = {._USE_STRIDED_LOAD = 0, ._USE_INDEXED_LOAD = 1, ._FINAL_N = 4, ._BARRETT_RED = 1, ._UNROLL_STOP = 3, ._FUSED_BUTTERFLY = 1};
    poly_mult_ntt_rvv_v3(dst, lhs, rhs, modulo, params_compressed);
}

/** Benchmark wrapper */
poly_mult_bench_result_t poly_mult_ntt_rvv_strided_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_strided, golden);
}

poly_mult_bench_result_t poly_mult_ntt_rvv_strided_barrett_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_strided_barrett, golden);
}

poly_mult_bench_result_t poly_mult_ntt_rvv_compressed_barrett_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_compressed_barrett, golden);
}

poly_mult_bench_result_t poly_mult_ntt_rvv_indexed_barrett_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_indexed_barrett, golden);
}

poly_mult_bench_result_t poly_mult_ntt_rvv_compressed_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_compressed, golden);
}

poly_mult_bench_result_t poly_mult_ntt_rvv_indexed_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_indexed, golden);
}

poly_mult_bench_result_t poly_mult_ntt_rvv_fastest_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, modulo, poly_mult_ntt_rvv_fastest, golden);
}