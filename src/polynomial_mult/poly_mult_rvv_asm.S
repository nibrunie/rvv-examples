// void rvv_ntt_transform_asm_helper(ntt_t* dst, int* coeffs, int _n, int level, int rootPowers[8][64]); 
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
// permutation
.loop_permute:
	vsetvli	a7, a3, e32, m8, ta, ma
	vle32.v	v8, (a5)
	vluxei32.v	v8, (a1), v8
	vse32.v	v8, (t0)
	sub	a3, a3, a7
	slli	a7, a7, 2
	add	a5, a5, a7
	add	t0, t0, a7
	bnez	a3, .loop_permute

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

// reconstruction level 6, 5, and 4
.loop_reconstruct_level_6_5_4:
	vsetvli	a7, a3, e32, m8, ta, mu
	// loading inputs into v8, vec_coeffs
	vle32.v	v8, (t0)
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

	li t1, 3329 // should be hoisted outside the loop
	vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction

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

	// assume t1 == 3329
	vrem.vx v8, v8, t1 // to be replaced by Barrett's reduction

	// storing results
	vse32.v	v8, (t0)
	sub	a3, a3, a7
	slli	a7, a7, 2
	add	t0, t0, a7
	// add t2, t2, a7 //  <= rootPowers address increment should not be done (each row only contains 64-element)
	// add t3, t3, a7 //  <= rootPowers address increment should not be done (each row only contains 64-element)
	bnez	a3, .loop_reconstruct_level_6_5_4
# %bb.37:                               # %._crit_edge111.3
	ld	t0, 8(a0) // t0 <- destination buffer
	ret
.Lfunc_end3:
	.size	rvv_ntt_transform_asm_internal, .Lfunc_end3-rvv_ntt_transform_asm_internal
                                        # -- End function