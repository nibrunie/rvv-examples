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

// building mask
	li t0, 0xaa
	vsetvli	a6, zero, e8, m1, ta, ma
	vmv.v.x v1, t0

	li t0, 0xcc
	vmv.v.x v2, t0

	li t0, 0xf0
	vmv.v.x v3, t0

	mv	a3, a2 // a3 <- avl
	ld	t0, 8(a0) // t0 <- destination buffer

	// int rootPowers[8][64]
	li t2, (5 * 64 * 4)
	add t2, a4, t2 // t2 <- rootPowers[5]

	li t3, (4 * 64 * 4)
	add t3, a4, t3 // t3 <- rootPowers[4]

	// hoisting constants used for modulo reduction (the second
	// one is only used when Barrett's method is used)
	li t1, 3329 // should be hoisted outside the loop
#ifndef USE_VREM_MODULO
	li t4, 5039
#endif

// reconstruction level 6, 5, and 4
.loop_reconstruct_level_6_5_4:
	vsetvli	a7, a3, e32, m8, ta, mu
	// loading indices
	vle32.v	v8, (a5)
	// performning permutation
	vluxei32.v	v8, (a1), v8
	// computing swapped elements
	vslidedown.vi v16, v8, 1
	vmv1r.v v0, v1 // v0 <- v1 (0xaa mask)
	vslideup.vi v16, v8, 1, v0.t
	vrsub.vi v8, v8, 0, v0.t // negate
	vadd.vv v8, v16, v8

	// butterfly
    // swapping odd/even pairs of coefficients
	vle32.v v24, (t2) // loading twiddleFactor
	vmv1r.v v0, v2
    // vec_coeffs = FUNC_LMUL_MASKED(__riscv_vmul_vv_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vec_twiddleFactor, vl);
	vmul.vv v8, v8, v24, v0.t
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
	vsetvli	a6, x0, e32, m4, ta, mu
	vwmul.vx v24, v8, t4
	vnsra.wi v24, v24, 24
	vnmsac.vx v8, t1, v24

	vwmul.vx v16, v12, t4
	vnsra.wi v16, v16, 24
	vnmsac.vx v12, t1, v16
	vsetvli	a6, x0, e32, m8, ta, mu
	vmsge.vx v0, v8, t1
	vsub.vx v8, v8, t1, v0.t
#endif

	// butterfly
    // swapping odd/even pairs of coefficients
	vle32.v v24, (t3) // loading twiddleFactor
	vmv1r.v v0, v3
    // vec_coeffs = FUNC_LMUL_MASKED(__riscv_vmul_vv_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vec_twiddleFactor, vl);
	vmul.vv v8, v8, v24, v0.t
    // TYPE_LMUL(vint32) vec_swapped_coeffs = FUNC_LMUL(__riscv_vslidedown_vx_i32)(vec_coeffs, n / 2, vl);
	vslidedown.vi v16, v8, 4
    // vec_swapped_coeffs = FUNC_LMUL_MASKED(__riscv_vslideup_vx_i32)(mask_up_b4, vec_swapped_coeffs, vec_coeffs, n / 2, vl);
	vslideup.vi v16, v8, 4, v0.t
    // vec_coeffs = FUNC_LMUL_MASKED(__riscv_vneg_v_i32)(mask_up_b4, vec_coeffs, vec_coeffs, vl);
	vrsub.vi v8, v8, 0, v0.t
    // vec_coeffs = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_coeffs, vec_swapped_coeffs, vl);
	vadd.vv v8, v16, v8

	// assume t1 == 3329 and t4 == 5039
#ifdef USE_VREM_MODULO
	vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction
#else
	vsetvli	a6, x0, e32, m4, ta, mu
	vwmul.vx v24, v8, t4
	vnsra.wi v24, v24, 24
	vnmsac.vx v8, t1, v24
	vmsge.vx v0, v8, t1
	vsub.vx v8, v8, t1, v0.t

	vwmul.vx v16, v12, t4
	vnsra.wi v16, v16, 24
	vnmsac.vx v12, t1, v16
	vmsge.vx v0, v12, t1
	vsub.vx v12, v12, t1, v0.t
	vsetvli	a6, x0, e32, m8, ta, mu
#endif

	// storing results
	vse32.v	v8, (t0)
	sub	a3, a3, a7
	slli	a7, a7, 2
	add	t0, t0, a7
	add	a5, a5, a7
	// add t2, t2, a7 //  <= rootPowers address increment should not be done (each row only contains 64-element)
	// add t3, t3, a7 //  <= rootPowers address increment should not be done (each row only contains 64-element)
	bnez	a3, .loop_reconstruct_level_6_5_4

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
	li t5, 3 // t5 <= local_level
	li t6, 16 // t6 <= n

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
    //         for (size_t vl; avl > 0; avl -= vl, even_coeffs += vl, odd_coeffs += vl, twiddleFactor += vl)
    //         {
.ntt_level_avl_loop:
    //             vl = FUNC_LMUL(__riscv_vsetvl_e32)(avl);
	vsetvli	a7, a3, e32, m4, ta, ma
    //             // TODO even coefficients should be loaded after odd+twiddleFactor as they are needed last
    //             TYPE_LMUL(vint32) vec_odd_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) odd_coeffs, vl);
	vle32.v	v16, (odd_coeffs_addr)
    //             TYPE_LMUL(vint32) vec_twiddleFactor = FUNC_LMUL(__riscv_vle32_v_i32)((int*) twiddleFactor, vl);
	vle32.v	v24, (t4)
    //             TYPE_LMUL(vint32) vec_odd_results = FUNC_LMUL(__riscv_vmul_vv_i32)(vec_odd_coeffs, vec_twiddleFactor, vl);
	vmul.vv v16, v16, v24
    //             TYPE_LMUL(vint32) vec_even_coeffs = FUNC_LMUL(__riscv_vle32_v_i32)((int*) even_coeffs, vl);
	vle32.v	v8, (even_coeffs_addr)
    //             TYPE_LMUL(vint32) vec_even_results = FUNC_LMUL(__riscv_vadd_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
	vadd.vv v24, v8, v16
    //             vec_odd_results = FUNC_LMUL(__riscv_vsub_vv_i32)(vec_even_coeffs, vec_odd_results, vl);
	vsub.vv v16, v8, v16

#ifdef USE_VREM_MODULO
	li t2, 3329
    //                 // even results
    //                 vec_even_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_even_results, dst->modulo, vl);
	vrem.vx v24, v24, t2 // to be replaced by Barrett's reduction
    //                 // odd results
    //                 vec_odd_results = FUNC_LMUL(__riscv_vrem_vx_i32)(vec_odd_results, dst->modulo, vl);
	vrem.vx v16, v16, t2 // to be replaced by Barrett's reduction
#else
    //                 vec_odd_results = rvv_barrett_reduction(vec_odd_results, vl);
    //                 vec_even_results = rvv_barrett_reduction(vec_even_results, vl);
	// Barrett's reduction of v24 / vec_even_results
	// TODO: maximum of LMUL=4 (because widening requires EMUL=8), could
	// be optmized with vl/LMUL change
	li t2, 5039
	vwmul.vx v8, v24, t2
	vnsra.wi v8, v8, 24
	li t2, 3329
	vnmsac.vx v24, t2, v8
	vmsge.vx v0, v24, t2
	vsub.vx v24, v24, t2, v0.t

	// Barrett's reduction of v16 / vec_odd_results
	li t2, 5039
	vwmul.vx v8, v16, t2
	vnsra.wi v8, v8, 24
	li t2, 3329
	vnmsac.vx v16, t2, v8
	vmsge.vx v0, v16, t2
	vsub.vx v16, v16, t2, v0.t
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
