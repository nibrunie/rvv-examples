// void rvv_ntt_transform_asm_internal(ntt_t* dst, int* coeffs, int _n, int level, int rootPowers[8][64]); 
// a0: destination buffer
// a1: input buffer
// a2: number of elements (should be 128)
// a3: level (should be )
// a4: address for root powers 2D array
.globl	rvv_ntt_transform_asm_internal # -- Begin function rvv_ntt_transform_asm_internal
	.p2align	1
	.type	rvv_ntt_transform_asm_internal,@function
rvv_ntt_transform_asm_internal:       # @rvv_ntt_transform_asm_internal
	ld	t0, 8(a0)
# %bb.1:                                # %.lr.ph.preheader
.Lpcrel_hi3:
	auipc	a3, %got_pcrel_hi(ntt_coeff_indices_128)
	ld	a5, %pcrel_lo(.Lpcrel_hi3)(a3)
	mv	a3, a2

// building mask and twiddle factors
    // v1: mask for level 6
	// v2: mask for level 5
	// v3: twiddle factors for level 4
	// v4v5: twiddle factors for level 3
	// v24v25v26v27: twiddle factors for level 2
	li t0, 0xaa
	vsetvli	a6, zero, e8, m1, ta, ma
	vmv.v.x v1, t0

	li t0, 0xcc
	vmv.v.x v2, t0

	mv	a3, a2 // a3 <- avl
	ld	t0, 8(a0) // t0 <- destination buffer

	// int rootPowers[8][64]
	li t2, (5 * 64 * 4)
	add t2, a4, t2 // t2 <- rootPowers[5]

	li t3, (4 * 64 * 4)
	add t3, a4, t3 // t3 <- rootPowers[4]
	vl1re32.v v3, (t3) // twiddle factors for level 4

	li t5, (3 * 64 * 4)
	add t5, a4, t5 // t5 <- rootPowers[3]
	vl2re32.v v4, (t5) // twiddle factors for level 3

	li t6, (2 * 64 * 4)
	add t6, a4, t6 // t6 <- rootPowers[2]
	vl4re32.v v24, (t6) // twiddle factors for level 2

	// hoisting constants used for modulo reduction (the second
	// one is only used when Barrett's method is used)
	li t1, 3329 // should be hoisted outside the loop
#ifndef USE_VREM_MODULO
	// t4 <- 1290167 = 2^32 / a(was 5039 = 2^24 / q)
	li t4, 1290167
#endif

// reconstruction level 6, 5, 4, 3 and 2
.loop_reconstruct_level_6_5_4_3_2:
	vsetvli	a7, a3, e32, m8, ta, mu
	// loading indices
	vle32.v	v8, (a5)
	// performing permutation
	vluxei32.v	v8, (a1), v8

	// level 6
	// computing swapped elements
	vslidedown.vi v16, v8, 1
	vmv1r.v v0, v1 // v0 <- v1 (0xaa mask)
	vslideup.vi v16, v8, 1, v0.t
	vrsub.vi v8, v8, 0, v0.t // negate
	vadd.vv v8, v16, v8

	// level 5 butterfly
    // swapping odd/even pairs of coefficients
	vle32.v v16, (t2) // loading twiddleFactor
	vmv1r.v v0, v2
    // vec_coeffs = FUNC_LMUL_MASKED(__riscv_vmul_vv_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vec_twiddleFactor, vl);
	vmul.vv v8, v8, v16, v0.t
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
	vmsge.vx v0, v8, t1
	vsub.vx v8, v8, t1, v0.t
#endif

	// level 4
	// butterfly
	// a single 4-element group of twiddle factor is required
	// here we assume VLEN=128 and perform the coefficient group swap directly
	// by using register indexing
	vsetvli	a6, x0, e32, m1, ta, ma
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


	vsetvli	a6, x0, e32, m8, ta, mu
	// assume t1 == 3329 and t4 == 1290167
#ifdef USE_VREM_MODULO
	vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
#else
	vmulh.vx v16, v8, t4
	vnmsac.vx v8, t1, v16
	vmsge.vx v0, v8, t1
	vsub.vx v8, v8, t1, v0.t
#endif


	// level 3
	// butterfly
	// a single 8-element group of twiddle factor is required
	vsetvli	a6, x0, e32, m2, ta, ma

	vmul.vv v18, v10, v4
	vmul.vv v22, v14, v4

	vsub.vv v10, v8, v18
	vadd.vv v8, v8, v18

	vsub.vv v14, v12, v22
	vadd.vv v12, v12, v22

	vsetvli	a6, x0, e32, m8, ta, mu
	// assume t1 == 3329 and t4 == 1290167
#ifdef USE_VREM_MODULO
	vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
#else
	vmulh.vx v16, v8, t4
	vnmsac.vx v8, t1, v16
	vmsge.vx v0, v8, t1
	vsub.vx v8, v8, t1, v0.t
#endif

	// level 2
	// butterfly
	// a single 16-element group of twiddle factor is required
	vsetvli	a6, x0, e32, m4, ta, ma

	vmul.vv v20, v12, v24

	vsub.vv v12, v8, v20
	vadd.vv v8, v8, v20

	// assume t1 == 3329 and t4 == 1290167
	vsetvli	a6, x0, e32, m8, ta, ma
#ifdef USE_VREM_MODULO
	vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
#else
	vmulh.vx v16, v8, t4
	vnmsac.vx v8, t1, v16
	vmsge.vx v0, v8, t1
	vsub.vx v8, v8, t1, v0.t
#endif

	// storing results
	vse32.v	v8, (t0)
	sub	a3, a3, a7
	slli	a7, a7, 2
	add	t0, t0, a7
	add	a5, a5, a7
	bnez	a3, .loop_reconstruct_level_6_5_4_3_2

// 
// last generic levels
#define even_coeffs_addr t3
#define odd_coeffs_addr t1
#define twiddle_factors_addr t4
//   t0: coeffs_a
//   t1: odd_coeffs
//   t2: temporary
//   t3: even_coeffs
//   t4: twiddleFactor
//   t5: local_level
//   t6: n
//   a5: j
//   a3: avl
//
//   half_n is not materialized (we use t6 >> 1 instead)
//   t2 is used as a temporary register whenever one is need
    //    n = 16;
    //    local_level = 3;
	li t5, 1 // t5 <= local_level
	li t6, 64 // t6 <= n

    // for (; local_level >= 0; n = 2 * n, local_level--) {
.ntt_level:
	// int j = 0
	li a5, 0

.ntt_j_loop:
    //     for (int j = 0; j < m; j++) {
    //         size_t avl = half_n;
	srli a3, t6, 1 // a3 <= avl = half_n
	// 
	// t0 <- coeffs_a
	ld	even_coeffs_addr, 8(a0)
	// 
    //         int* even_coeffs = coeffs_a + 2 * j * half_n;
	mul a7, a5, t6 // 2 * j * half_n = j * n
	sll a7, a7, 2 // sizeof(int)=4 * 2 * j * half_n
	add even_coeffs_addr, even_coeffs_addr, a7 // t0 <- even_coeffs
	// 
    //         int* odd_coeffs = even_coeffs + half_n;
	slli a7, t6, 1 // sizeof(int) * half_n = 4 * n / 2 = 2 *n
	add odd_coeffs_addr, even_coeffs_addr, a7 // t1 <- odd_coeffs
    //         int* twiddleFactor = rootPowers[local_level];
	li twiddle_factors_addr, (64 * 4)
	mul twiddle_factors_addr, twiddle_factors_addr, t5 
	add twiddle_factors_addr, a4, twiddle_factors_addr // t4 <- rootPowers[local_level]
	//
	// q (used for modulo reduction)
	li a6, 3329
	// t2 <- 1290167 = 2^32 / a(was 5039 = 2^24 / q)
	li t2, 1290167
    //         for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, twiddleFactor += vl)
    //         {
.ntt_level_avl_loop:
    //             vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
	vsetvli	a7, a3, e32, m8, ta, ma
    //             // TODO even coefficients should be loaded after odd+twiddleFactor as they are needed last
    //             TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) odd_coeffs, vl);
	vle32.v	v16, (odd_coeffs_addr)
    //             TYPE_LMUL(vint32) vec_twiddleFactor = FUNC_LMUL(__riscv_vle32_v_i32)((int*) twiddleFactor, vl);
	vle32.v	v24, (t4)
    //             TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) even_coeffs, vl);
	vle32.v	v8, (even_coeffs_addr)
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
    //                 // odd results
    //                 vec_odd_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_odd_results, dst->modulo, vl);
	vrem.vx v16, v16, a6 // to be replaced by Barrett's reduction
#else
    //                 vec_odd_results = rvv_barrett_reduction(vec_odd_results, vl);
    //                 vec_even_results = rvv_barrett_reduction(vec_even_results, vl);
	// Barrett's reduction of v24 / vec_even_results
	vmulh.vx v8, v24, t2
	vnmsac.vx v24, a6, v8
	vmsge.vx v0, v24, a6
	vsub.vx v24, v24, a6, v0.t

	// Barrett's reduction of v16 / vec_odd_results
	vmulh.vx v8, v16, t2
	vnmsac.vx v16, a6, v8
	vmsge.vx v0, v16, a6
	vsub.vx v16, v16, a6, v0.t
#endif
    //             }
    //             FUNC_LMUL(__riscv_vse32_v_i32)(even_coeffs, vec_even_results, vl);
	vse32.v	v24, (even_coeffs_addr)
    //             FUNC_LMUL(__riscv_vse32_v_i32)(odd_coeffs, vec_odd_results, vl);
	vse32.v	v16, (odd_coeffs_addr)

	sub	a3, a3, a7
	slli	a7, a7, 2
	add	even_coeffs_addr, even_coeffs_addr, a7
	add odd_coeffs_addr, odd_coeffs_addr, a7
	add t4, t4, a7 // twiddle_factors += vl
	bnez	a3, .ntt_level_avl_loop
    //         }
	addi a5, a5, 1 // j++
    //     const int m = 1 << local_level;
	li t2, 1
	sll t2, t2, t5 // t2 <= m = (1 << local_level)
	sub a6, a5, t2 // j - m
	bnez a6, .ntt_j_loop
    //     } 
    // }
	addi t5, t5, -1 // local_level --
	sll t6, t6, 1 // n = 2 * n
	bgez t5, .ntt_level

	ret
.Lfunc_end3:
	.size	rvv_ntt_transform_asm_internal, .Lfunc_end3-rvv_ntt_transform_asm_internal
                                        # -- End function

// void rvv_ntt_mult_scale_asm(int* dst, int* lhs, int* rhs); 
// a0: destination buffer
// a1: lhs input buffer
// a2: rhs input buffer
.globl	rvv_ntt_mult_scale_asm # -- Begin function rvv_ntt_transform_asm_internal
	.p2align	1
	.type	rvv_ntt_mult_scale_asm,@function
rvv_ntt_mult_scale_asm:       # @rvv_ntt_transform_asm_internal
	li a3, 128 // a3 = avl = 128
	// q (used for modulo reduction)
	li a6, 3329
	// t2 <- 1290167 = 2^32 / a(was 5039 = 2^24 / q)
	li t2, 1290167
	// 1 / n mod q
	li t3, 3303
.rvv_ntt_mult_scale_avl_loop:
	vsetvli	a7, a3, e32, m8, ta, ma
	// loading lhs coefficients
	vle32.v	v16, (a1)
	// loading rhs coefficients
	vle32.v	v8, (a2)

	// element-wise multiplication
	vmul.vv v16, v16, v8

	// modulo reduction for element wise multiplication
#ifdef USE_VREM_MODULO
	vrem.vx v16, v16, a6 // to be replaced by Barrett's reduction
#else
	// Barrett's reduction of v16 vector register group
	vmulh.vx v8, v16, t2
	vnmsac.vx v16, a6, v8
	vmsge.vx v0, v16, a6
	vsub.vx v16, v16, a6, v0.t
#endif

	// degree scaling
	vmul.vx v16, v16, t3

	// modulo reduction for element wise multiplication
#ifdef USE_VREM_MODULO
	vrem.vx v16, v16, a6 // to be replaced by Barrett's reduction
#else
	// Barrett's reduction of v16 vector register group
	vmulh.vx v8, v16, t2
	vnmsac.vx v16, a6, v8
	vmsge.vx v0, v16, a6
	vsub.vx v16, v16, a6, v0.t
#endif

	vse32.v	v16, (a0)

	sub	a3, a3, a7
	slli	a7, a7, 2
	add a0, a0, a7 // updating destination address
	add	a1, a1, a7 // updating lhs address
	add a2, a2, a7 // updating rhs address
	bnez	a3, .rvv_ntt_mult_scale_avl_loop
.rvv_ntt_mult_scale_asm_end:
	ret
	.size	rvv_ntt_mult_scale_asm, .rvv_ntt_mult_scale_asm_end-rvv_ntt_transform_asm_internal