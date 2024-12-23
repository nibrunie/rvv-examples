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
	vsetvli	a6, zero, e8, m8, ta, ma
	vmv.v.x v0, t0

	mv	a3, a2 // a3 <- avl
	ld	t0, 8(a0) // t0 <- destination buffer

// reconstruction level 6, 5, and 4
# %bb.37:                               # %._crit_edge111.3
	ret
.Lfunc_end3:
	.size	rvv_ntt_transform_asm_internal, .Lfunc_end3-rvv_ntt_transform_asm_internal
                                        # -- End function