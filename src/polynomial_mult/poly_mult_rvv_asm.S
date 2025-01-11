// void rvv_ntt_transform_asm_internal(int* dst, int* coeffs, int rootPowers[8][64], int finalCorrection); 
// a0: destination buffer
// a1: input buffer
// a2: address for root powers 2D array
// a3: perform barrett's modulo reduction final correction (in level 0)
.globl    rvv_ntt_transform_asm_internal # -- Begin function rvv_ntt_transform_asm_internal
    .p2align    1
    .type    rvv_ntt_transform_asm_internal,@function
rvv_ntt_transform_asm_internal:       # @rvv_ntt_transform_asm_internal
// materializing the permutation coefficient indices array address (in register a5)
.Lpcrel_hi3:
    auipc    a4, %got_pcrel_hi(ntt_coeff_indices_128)
    ld    a5, %pcrel_lo(.Lpcrel_hi3)(a4)

// materializing the consntant application vector length (in register a4)
    li    a4, 128    // a4 = avl = 128

// macros DISABLE_* are used to selectively disable some code sections.
// This makes the code non functional but can be used to evaluate the
// latency impacts of different selections of sections.
#ifndef DISABLE_NTT_FIRST_FUSED_LEVELS
#ifndef DISABLE_NTT_FUSED_LEVELS_INIT

// building mask and twiddle factors for level 6 through 2
// static allocation of vector registers (invariant in loop)
//   v1: mask for level 6
//   v2: mask for level 5
//   v3: twiddle factors for level 4
//   v4v5v6v7: twiddle factors for level 2
//   v24v25v26v27v28v29v30v31: twiddle factors for level 5

// mask for level 6 (v1)
    li t0, 0xaa
    vsetvli    a6, zero, e8, m1, ta, ma
    vmv.v.x v1, t0

// mask for level 5 (v2)
    li t0, 0xcc
    vmv.v.x v2, t0

    mv    t0, a0 // t0 <- destination buffer

// increment update (in bytes) for array addresses (in t6)
    vsetvli    a7, a4, e32, m8, ta, mu
    slli t6, a7, 2 // sizeof(int) * vl => byte address increment update

// computing addresses for rootPowers (load from int rootPowers[8][64])
// + pre-loading some vector register groups

// twiddle factors for level 2 (in v4v5v6v7)
// The use of a whole vector register group vector load prevents us from
// needing to insert a vsetvl* instruction here.
    li t5, (2 * 64 * 4)
    add t5, a2, t5 // t6 <- rootPowers[2]
    vl4re32.v v4, (t5) // twiddle factors for level 2

// twiddle factors for level 4 (in v3)
    li t3, (4 * 64 * 4)
    add t3, a2, t3 // t3 <- rootPowers[4]
    vl1re32.v v3, (t3) // twiddle factors for level 4

// twiddle factors for level 5 (in v24v25v26v27v28v29v30v31)
// Since we used a masked load we cannot rely on a whole vector register group load
// and have to make sure vtype is properly configured (LMUL=8, VLMAX)
    li t2, (5 * 64 * 4)
    add t2, a2, t2 // t2 <- rootPowers[5]
    vmv.v.i v24, 1
    vmv1r.v v0, v2
    vle32.v v24, (t2), v0.t // loading twiddleFactor for level 5

// last address computation since we need t5 to keep the address of rootPowers[2]
// t5 contains the address of the twiddle factors array for level 3
    li t5, (3 * 64 * 4)
    add t5, a2, t5 // t5 <- rootPowers[3]

// hoisting constants used for modulo reduction (the second
// one is only used when Barrett's method is used)
    li t1, 3329 // should be hoisted outside the loop


// The USE_VREM_MODULO macro is used to switch between the use of the
// vrem instruction and the use of the Barrett's reduction method.
// Barrett's reduction is always faster in our experiments but vrem
// is useful to debug the implementation.
#ifndef USE_VREM_MODULO
    // t4 <- 1290167 = 2^32 / q (was 5039 = 2^24 / q)
    li t4, 1290167
#endif // USE_VREM_MODULO

#else // DISABLE_NTT_FUSED_LEVELS_INIT
    li    t6, 128 // a4 <- avl
    mv    t0, a0 // t0 <- destination buffer
#endif // DISABLE_NTT_FUSED_LEVELS_INIT

// reconstruction for levels 6, 5, 4, 3 and 2
// register mapping in this loop
//   a1: start index of input coefficient buffer
//   a4: avl
//   a5: current address of index buffer
//   a7: vl
//
//   t0: current address of destination buffer
//   t1: q
//   t2: current address of twiddle factor for level 5
//   t4: 2^32 // q
//   t5: address of twiddle factor for level 3 (loop invariant)
//   t6: sizeof(int) * vl => byte address increment update
.loop_reconstruct_level_6_5_4_3_2:
#ifndef DISABLE_NTT_FUSED_LEVELS_PERMUTE
    // loading indices
    vle32.v    v8, (a5)
    // performing permutation
    vluxei32.v    v8, (a1), v8

    vmv1r.v v0, v1 // v0 <- v1 (0xaa mask)

#endif // DISABLE_NTT_FUSED_LEVELS_PERMUTE
#ifndef DISABLE_NTT_FUSED_LEVELS_LVL6
    // level 6
    // computing swapped elements
    vslidedown.vi v16, v8, 1
    vslideup.vi v16, v8, 1, v0.t
    vrsub.vi v8, v8, 0, v0.t // negate
    vadd.vv v8, v16, v8
#endif // DISABLE_NTT_FUSED_LEVELS_LVL6

#ifndef DISABLE_NTT_FUSED_LEVELS_LVL5
    // level 5 butterfly
    // swapping odd/even pairs of coefficients
    // vec_coeffs = FUNC_LMUL_MASKED(__riscv_vmul_vv_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vec_twiddleFactor, vl);
    vmul.vv v8, v8, v24
    vmv1r.v v0, v2
    // TYPE_LMUL(vint32) vec_swapped_coeffs = FUNC_LMUL(__riscv_vslidedown_vx_i32)(vec_coeffs, n / 2, vl);
    vslidedown.vi v16, v8, 2
    // vec_swapped_coeffs = FUNC_LMUL_MASKED(__riscv_vslideup_vx_i32)(mask_up_b4, vec_swapped_coeffs, vec_coeffs, n / 2, vl);
    vslideup.vi v16, v8, 2, v0.t
    // vec_coeffs = FUNC_LMUL_MASKED(__riscv_vneg_v_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vl);
    vrsub.vi v8, v8, 0, v0.t
    // vec_coeffs = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_coeffs, vec_swapped_coeffs, vl);
    vadd.vv v8, v16, v8

#ifdef USE_VREM_MODULO
    vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
#else
    vmulh.vx v16, v8, t4
    vnmsac.vx v8, t1, v16
#endif

#endif // DISABLE_NTT_FUSED_LEVELS_LVL5

#ifndef DISABLE_NTT_FUSED_LEVELS_LVL4
    // level 4
    // butterfly
    // a single 4-element group of twiddle factor is required
    // here we assume VLEN=128 and perform the coefficient group swap directly
    // by using register indexing
    vsetvli    a6, x0, e32, m1, ta, ma
    // explicit 4-element swap assuming VLEN=128
    vmul.vv v17, v9, v3
    vmul.vv v19, v11, v3
    vmul.vv v21, v13, v3
    vmul.vv v23, v15, v3

    vsub.vv v9, v8, v17
    vadd.vv v8, v8, v17

    vsub.vv v11, v10, v19
    vadd.vv v10, v10, v19

    vsub.vv v13, v12, v21
    vadd.vv v12, v12, v21

    vsub.vv v15, v14, v23
    vadd.vv v14, v14, v23


    vsetvli    a6, x0, e32, m8, ta, ma
    // assume t1 == 3329 and t4 == 1290167
    // loading of twiddle factors for level 3 in v16v17 is anticipated
#ifdef USE_VREM_MODULO
    vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
    vl2re32.v v16, (t5) // twiddle factors for level 3
#else
    vmulh.vx v16, v8, t4
    vnmsac.vx v8, t1, v16
    vl2re32.v v16, (t5) // twiddle factors for level 3
#endif

#endif // DISABLE_NTT_FUSED_LEVELS_LVL4

#ifndef DISABLE_NTT_FUSED_LEVELS_LVL3
    // level 3 butterfly
    // a single 8-element group of twiddle factor is required (assumed to be present in v4)
    vsetvli    a6, x0, e32, m2, ta, ma

    vmul.vv v18, v10, v16
    vmul.vv v22, v14, v16

    vsub.vv v10, v8, v18
    vadd.vv v8, v8, v18

    vsub.vv v14, v12, v22
    vadd.vv v12, v12, v22

    vsetvli    a6, x0, e32, m8, ta, ma
    // assume t1 == 3329 (q) and t4 == 1290167 (2^32 // q)
#ifdef USE_VREM_MODULO
    vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
#else
    vmulh.vx v16, v8, t4
    vnmsac.vx v8, t1, v16
#endif

#endif // DISABLE_NTT_FUSED_LEVELS_LVL3

#ifndef DISABLE_NTT_FUSED_LEVELS_LVL2
    // level 2 // butterfly
    // a single 16-element group of twiddle factor is required (assumed to be present in v24)
    vsetvli    a6, x0, e32, m4, ta, ma

    vmul.vv v20, v12, v4

    vsub.vv v12, v8, v20
    vadd.vv v8, v8, v20

    // assume t1 == 3329 (q) and t4 == 1290167 (2^32 // q)
    vsetvli    a7, a4, e32, m8, ta, ma
#ifdef USE_VREM_MODULO
    vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
#else
    vmulh.vx v16, v8, t4
    vnmsac.vx v8, t1, v16
#endif

#endif // DISABLE_NTT_FUSED_LEVELS_LVL2

#ifndef DISABLE_NTT_FUSED_LEVELS_STORING
    // storing results
    vse32.v    v8, (t0)

#endif // DISABLE_NTT_FUSED_LEVELS_STORING

    // updating addresses and loop counter
    addi    a4, a4, -32
    add     a5, a5, t6
    add     t0, t0, t6
    bnez    a4, .loop_reconstruct_level_6_5_4_3_2

#endif // DISABLE_NTT_FIRST_FUSED_LEVELS

// 
// last generic levels
#define even_coeffs_addr t3
#define odd_coeffs_addr t1
#define twiddle_factors_addr t4
//   t1: odd_coeffs address address
//   t2: q(3329) / temporary outside innermost loop
//   t3: even_coeffs current address
//   t4: twiddleFactor current address
//   t5: local_level
//   t6: n
//   a4: avl
//   a5: j
//   a6: 2^32 // Q
//
//   half_n is not materialized (we use t6 >> 1 instead)
//   t2 is used as a temporary register whenever one is need
    //    n = 16;
    //    local_level = 3;
    li t5, 1 // t5 <= local_level
    li t6, 64 // t6 <= n
    // q (used for modulo reduction)
    li a6, 3329
    //
    // a7 <- 1290167 = 2^32 / q (was 5039 = 2^24 / q)
    li a7, 1290167

    // constant for address update
    li t0, 128 // sizeof(int) * vl => byte address increment update
    // for (; local_level >= 0; n = 2 * n, local_level--) {
#ifndef DISABLE_NTT_LEVEL_1
.ntt_level:
    // int j = 0
    li a5, 0

    li a4, 32 // minimum avl to get vl = vlmax =32 (for e32 and VLEN=128)
    // setting outside the loop since vl should be uniform for all iterations
    // mask policy for level 0 must be mu (mask-undisturbed) so we need to set mu here
    vsetvli    x0, a4, e32, m8, ta, mu
// level 1
.ntt_j_loop:
    //     for (int j = 0; j < m; j++) {
    //         size_t avl = half_n;
    li a4, 32 // a4 <= avl = half_n = t6 / 2 = 32
    // 
    // t3 <- coeffs_a
    mv    even_coeffs_addr, a0
    // 
    //         int* even_coeffs = coeffs_a + 2 * j * half_n;
    mul t2, a5, t6 // 2 * j * half_n = j * n
    // sll a7, a7, 2 // sizeof(int)=4 * 2 * j * half_n
    sh2add even_coeffs_addr, t2, even_coeffs_addr // t3 <- even_coeffs
    // 
    //         int* odd_coeffs = even_coeffs + half_n;
    sh1add odd_coeffs_addr, t6, even_coeffs_addr // t1 <- odd_coeffs = evens_coeffs + sizeof(int) * half_n = 4 * n / 2 = 2 *n
    //         int* twiddleFactor = rootPowers[local_level];
    // rootPowers offset for local_level=1 is 256
    addi twiddle_factors_addr, a2, 256 // t4 <- rootPowers[local_level=1]
    //
    //         for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, twiddleFactor += vl)
    //         {
    // There are only 32 twiddle factors for level 1 so a single LMUL=8 vle32 is enough
    // and can be factorized across the two AVL loop
    //             TYPE_LMUL(vint32) vec_twiddleFactor = FUNC_LMUL(__riscv_vle32_v_i32)((int*) twiddleFactor, vl);
    vle32.v    v0, (t4)
    //             vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
.ntt_level_avl_loop:
    //             TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) odd_coeffs, vl);
    vle32.v    v16, (odd_coeffs_addr)
    //             TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) even_coeffs, vl);
    vle32.v    v8, (even_coeffs_addr)
    //             TYPE_LMUL(vint32) vec_odd_results = FUNC_LMUL(__riscv_vmul_vv_i32)(vec_odd_coeffs, vec_twiddleFactor, vl);
    vmul.vv v16, v16, v0
    //             TYPE_LMUL(vint32) vec_even_results = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
    vadd.vv v24, v8, v16
    //             vec_odd_results = FUNC_LMUL(__riscv_vsub_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
    vsub.vv v16, v8, v16

#ifdef USE_VREM_MODULO
    //                 // even results
    //                 vec_even_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_even_results, dst->modulo, vl);
    vrem.vx v24, v24, a6 // to be replaced by Barrett's reduction
    //             }
    //             FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);
    vse32.v    v24, (even_coeffs_addr)
    //                 // odd results
    //                 vec_odd_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_odd_results, dst->modulo, vl);
    vrem.vx v16, v16, a6 // to be replaced by Barrett's reduction
#else
    //                 vec_odd_results = rvv_barrett_reduction(vec_odd_results, vl);
    //                 vec_even_results = rvv_barrett_reduction(vec_even_results, vl);
    // Barrett's reduction of v24 / vec_even_results
    vmulh.vx v8, v24, a7
    vnmsac.vx v24, a6, v8
    //             }
    //             FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);
    vse32.v    v24, (even_coeffs_addr)

    // Barrett's reduction of v16 / vec_odd_results
    vmulh.vx v8, v16, a7
    vnmsac.vx v16, a6, v8
#endif
    //             FUNC_LMUL(__riscv_vse32_v_i32)(odd_coeffs, vec_odd_results, vl);
    vse32.v    v16, (odd_coeffs_addr)

    addi    a4, a4, -32
    add     even_coeffs_addr, even_coeffs_addr, t0
    add     odd_coeffs_addr, odd_coeffs_addr, t0
    bnez    a4, .ntt_level_avl_loop
    //         }
    addi a5, a5, 1 // j++
    //     const int m = 1 << local_level;
    li t2, 1
    sll t2, t2, t5 // t2 <= m = (1 << local_level)
    sub t2, a5, t2 // j - m
    bnez t2, .ntt_j_loop

#endif // DISABLE_NTT_LEVEL_1

#ifndef DISABLE_NTT_LEVEL_0
    //     } 
    // }
    sll t6, t6, 1 // n = 2 * n

    // level 0
    li a5, 0 // j <- 0
    bnez a3, .ntt_correction_level_0_j_loop
.ntt_level_0_j_loop:
    //     for (int j = 0; j < m; j++) {
    //         size_t avl = half_n;
    li a4, 64 // a4 <= avl = half_n = 64
    // 
    //         int* even_coeffs = coeffs_a + 2 * j * half_n;
    // t3 <- coeffs_a = even_coeffs_addr
    mv    even_coeffs_addr, a0
    // 
    //         int* odd_coeffs = even_coeffs + half_n;
    sh1add odd_coeffs_addr, t6, even_coeffs_addr // t1 <- odd_coeffs = evens_coeffs + sizeof(int) * half_n = 4 * n / 2 = 2 *n
    //         int* twiddleFactor = rootPowers[local_level];
    mv twiddle_factors_addr, a2 // t4 <- rootPowers[local_level=0]
    //
    //         for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, twiddleFactor += vl)
    //         {
    //             vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
    // vsetvli    x0, a4, e32, m8, ta, ma // setting outside the loop since vl should be uniform for all iterations
.ntt_level_0_avl_loop:
    //             TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) odd_coeffs, vl);
    vle32.v    v16, (odd_coeffs_addr)
    //             TYPE_LMUL(vint32) vec_twiddleFactor = FUNC_LMUL(__riscv_vle32_v_i32)((int*) twiddleFactor, vl);
    vle32.v    v24, (t4)
    //             TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) even_coeffs, vl);
    vle32.v    v8, (even_coeffs_addr)
    //             TYPE_LMUL(vint32) vec_odd_results = FUNC_LMUL(__riscv_vmul_vv_i32)(vec_odd_coeffs, vec_twiddleFactor, vl);
    vmul.vv v16, v16, v24
    //             TYPE_LMUL(vint32) vec_even_results = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
    vadd.vv v24, v8, v16
    //             vec_odd_results = FUNC_LMUL(__riscv_vsub_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
    vsub.vv v16, v8, v16

#ifdef USE_VREM_MODULO
    //                 // even results
    //                 vec_even_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_even_results, dst->modulo, vl);
    vrem.vx v24, v24, a6 // to be replaced by Barrett's reduction
    //             }
    //             FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);
    vse32.v    v24, (even_coeffs_addr)
    //                 // odd results
    //                 vec_odd_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_odd_results, dst->modulo, vl);
    vrem.vx v16, v16, a6 // to be replaced by Barrett's reduction
#else
    //                 vec_odd_results = rvv_barrett_reduction(vec_odd_results, vl);
    //                 vec_even_results = rvv_barrett_reduction(vec_even_results, vl);
    // Barrett's reduction of v24 / vec_even_results
    vmulh.vx v8, v24, a7
    vnmsac.vx v24, a6, v8
    //             }
    //             FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);
    vse32.v    v24, (even_coeffs_addr)

    // Barrett's reduction of v16 / vec_odd_results
    vmulh.vx v8, v16, a7
    vnmsac.vx v16, a6, v8
#endif
    //             FUNC_LMUL(__riscv_vse32_v_i32)(odd_coeffs, vec_odd_results, vl);
    vse32.v    v16, (odd_coeffs_addr)

    addi    a4, a4, -32
    add     even_coeffs_addr, even_coeffs_addr, t0
    add     odd_coeffs_addr, odd_coeffs_addr, t0
    add     t4, t4, t0 // twiddle_factors += vl
    bnez    a4, .ntt_level_0_avl_loop
    //         }
    //     } 
    // }

    ret
.ntt_correction_level_0_j_loop:
    li a4, 64 // a4 <= avl = half_n = 64
    mv    even_coeffs_addr, a0
    sh1add odd_coeffs_addr, t6, even_coeffs_addr // t1 <- odd_coeffs = evens_coeffs + sizeof(int) * half_n = 4 * n / 2 = 2 *n
    mv twiddle_factors_addr, a2 // t4 <- rootPowers[local_level=0]

.ntt_correction_level_0_avl_loop:
    vle32.v    v16, (odd_coeffs_addr)
    vle32.v    v24, (t4)
    vle32.v    v8, (even_coeffs_addr)

    vmul.vv v16, v16, v24
    vadd.vv v24, v8, v16
    vsub.vv v16, v8, v16

#ifdef USE_VREM_MODULO
    vrem.vx v24, v24, a6 // to be replaced by Barrett's reduction
    vse32.v    v24, (even_coeffs_addr)
    vrem.vx v16, v16, a6 // to be replaced by Barrett's reduction
#else
    vmulh.vx v8, v24, a7
    vnmsac.vx v24, a6, v8
    vmsge.vx v0, v24, a6
    vsub.vx v24, v24, a6, v0.t

    vse32.v    v24, (even_coeffs_addr)

    // Barrett's reduction of v16 / vec_odd_results
    vmulh.vx v8, v16, a7
    vnmsac.vx v16, a6, v8
    vmsge.vx v0, v16, a6
    vsub.vx v16, v16, a6, v0.t
#endif
    vse32.v    v16, (odd_coeffs_addr)

    addi    a4, a4, -32
    add     even_coeffs_addr, even_coeffs_addr, t0
    add     odd_coeffs_addr, odd_coeffs_addr, t0
    add     t4, t4, t0 // twiddle_factors += vl
    bnez    a4, .ntt_correction_level_0_avl_loop
#endif // DISABLE_NTT_LEVEL_0

    ret
.Lfunc_end3:
    .size    rvv_ntt_transform_asm_internal, .Lfunc_end3-rvv_ntt_transform_asm_internal
    # -- End function

// void rvv_ntt_mult_scale_asm(int* dst, int* lhs, int* rhs); 
// a0: destination buffer
// a1: lhs input buffer
// a2: rhs input buffer
// Note: This function does not perform correction at the end of the Barrett's modulo reduction
//       which means that result elements are in the [0; 2*q[ range instead of [0; q[
//       (which would have required extra instuctions and will be carried out by
//       the final inverse NTT anyway)
.globl    rvv_ntt_mult_scale_asm # -- Begin function rvv_ntt_transform_asm_internal
    .p2align    1
    .type    rvv_ntt_mult_scale_asm,@function
rvv_ntt_mult_scale_asm:       # @rvv_ntt_transform_asm_internal
#ifndef DISABLE_NTT_MULT_SCALE
    li a4, 128 // a4 = avl = 128
    // q (used for modulo reduction)
    li a6, 3329
    // t2 <- 1290167 = 2^32 / q (was 5039 = 2^24 / q)
    li t2, 1290167
    // 1 / n mod q (for degree scaling)
    li t3, 3303
    // (1 / n % q) * 2^32 / q
    li t4, 4261422943

    vsetvli    a7, a4, e32, m8, ta, ma // setting outside the loop since vl should be uniform for all iterations
    li t0, 128 // sizeof(int) * vl => byte address increment update
.rvv_ntt_mult_scale_avl_loop:
    // loading lhs coefficients
    vle32.v    v16, (a1)
    // loading rhs coefficients
    vle32.v    v8, (a2)

    // element-wise multiplication
    vmul.vv v16, v16, v8

    // modulo reduction for element wise multiplication
#ifdef USE_VREM_MODULO
    vrem.vx v16, v16, a6 // to be replaced by Barrett's reduction
#else
    // Barrett's reduction of v16 vector register group
    vmulh.vx v8, v16, t2
    vnmsac.vx v16, a6, v8
#endif

    // degree scaling
    vmul.vx v24, v16, t3

    // modulo reduction for degree scaling
    // first multiplication uses pre-computed b*2^32 // q = (3303 << 32) / 3329
    // as in degree scaling, this factor is constant. This does not reduce
    // the number of instructions but allows to use an instruction which can be
    // parallelized with the raw degree scaling multiplication above.
#ifdef USE_VREM_MODULO
    vrem.vx v16, v24, a6 // to be replaced by Barrett's reduction
#else
    // Barrett's reduction of v16 vector register group
    vmulhu.vx v8, v16, t4
    vnmsac.vx v24, a6, v8
#endif

    vse32.v    v24, (a0)

    addi    a4, a4, -32
    add     a0, a0, t0 // updating destination address
    add     a1, a1, t0 // updating lhs address
    add     a2, a2, t0 // updating rhs address
    bnez    a4, .rvv_ntt_mult_scale_avl_loop
#endif // DISABLE_NTT_MULT_SCALE
.rvv_ntt_mult_scale_asm_end:
    ret
    .size    rvv_ntt_mult_scale_asm, .rvv_ntt_mult_scale_asm_end-rvv_ntt_transform_asm_internal
