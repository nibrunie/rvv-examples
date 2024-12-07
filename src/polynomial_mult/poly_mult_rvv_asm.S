
	.globl	rvv_ntt_transform_asm_helper # -- Begin function rvv_ntt_transform_asm_helper
	.p2align	1
	.type	rvv_ntt_transform_asm_helper,@function
rvv_ntt_transform_asm_helper:       # @rvv_ntt_transform_asm_helper
# %bb.0:
	addi	sp, sp, -32
	sd	ra, 24(sp)                      # 8-byte Folded Spill
	sd	s0, 16(sp)                      # 8-byte Folded Spill
	sd	s1, 8(sp)                       # 8-byte Folded Spill
	li	a3, 1
	bge	a3, a2, .LBB3_16
# %bb.1:
	vsetvli	t0, zero, e32, m1, ta, ma
	li	a3, -86
	vsetvli	zero, t0, e8, m1, ta, ma
	ld	a5, 8(a0)
	vmv.v.x	v8, a3
	lui	a3, %hi(ntt_coeff_indices_128)
	addi	s0, a3, %lo(ntt_coeff_indices_128)
	mv	a3, a2
.LBB3_2:                                # =>This Inner Loop Header: Depth=1
	vsetvli	a6, a3, e32, m4, ta, ma
	vle32.v	v12, (s0)
	vluxei32.v	v12, (a1), v12
	vse32.v	v12, (a5)
	sub	a3, a3, a6
	slli	a6, a6, 2
	add	s0, s0, a6
	add	a5, a5, a6
	bnez	a3, .LBB3_2
# %bb.3:
	ld	a7, 8(a0)
	vsetvli	a6, zero, e32, m4, ta, ma
	vsetvli	zero, t0, e8, m1, ta, ma
	vmv.v.i	v9, -16
	addi	a0, a4, 1024
	vsetvli	a1, zero, e32, m4, ta, ma
	vle32.v	v12, (a0)
	li	a0, -52
	vsetvli	zero, t0, e8, m1, ta, ma
	vmv.v.x	v10, a0
	addi	a0, a4, 1280
	vsetvli	a1, zero, e32, m4, ta, ma
	vle32.v	v16, (a0)
	lui	a1, 1
	addi	a0, a1, 943
	addi	a3, a1, -767
	lui	a1, 1048575
	addi	a5, a1, 767
	mv	a1, a7
.LBB3_4:                                # =>This Inner Loop Header: Depth=1
	vsetvli	s0, a2, e32, m4, ta, mu
	vle32.v	v20, (a1)
	vslidedown.vi	v24, v20, 1
	vmv1r.v	v0, v8
	vslideup.vi	v24, v20, 1, v0.t
	vrsub.vi	v20, v20, 0, v0.t
	vadd.vv	v20, v20, v24
	vmv1r.v	v0, v10
	vmul.vv	v20, v20, v16, v0.t
	vslidedown.vi	v24, v20, 2
	vslideup.vi	v24, v20, 2, v0.t
	vrsub.vi	v20, v20, 0, v0.t
	vadd.vv	v24, v20, v24
	vwmul.vx	v0, v24, a0
	vnsra.wi	v20, v0, 24
	vnmsub.vx	v20, a3, v24
	vmslt.vx	v11, v20, a3
	vmnot.m	v0, v11
	vadd.vx	v20, v20, a5, v0.t
	vmv1r.v	v0, v9
	vmul.vv	v20, v20, v12, v0.t
	vslidedown.vi	v24, v20, 4
	vslideup.vi	v24, v20, 4, v0.t
	vrsub.vi	v20, v20, 0, v0.t
	vadd.vv	v20, v20, v24
	vwmul.vx	v24, v20, a0
	vnsra.wi	v4, v24, 24
	vnmsub.vx	v4, a3, v20
	vmslt.vx	v11, v4, a3
	vmnot.m	v0, v11
	vadd.vx	v4, v4, a5, v0.t
	vse32.v	v4, (a1)
	sub	a2, a2, s0
	slli	s0, s0, 2
	add	a1, a1, s0
	bnez	a2, .LBB3_4
# %bb.5:
	li	a0, 3
	bgeu	a0, a6, .LBB3_17
# %bb.6:
	li	a0, 8
	bltu	a6, a0, .LBB3_18
# %bb.7:                                # %.lr.ph.split.preheader.preheader
	li	t0, 16
	li	a0, 3
	lui	a1, 1
	addi	t6, a1, 943
	addi	a3, a1, -767
	lui	a1, 1048575
	addi	a5, a1, 767
.LBB3_8:                                # %.lr.ph.split.preheader
                                        # =>This Loop Header: Depth=1
                                        #     Child Loop BB3_11 Depth 2
                                        #       Child Loop BB3_12 Depth 3
	mv	a6, a0
	li	a0, 1
	sllw	t1, a0, a6
	bgtz	t1, .LBB3_10
# %bb.9:                                # %.lr.ph.split.preheader
                                        #   in Loop: Header=BB3_8 Depth=1
	li	t1, 1
.LBB3_10:                               # %.lr.ph.split.preheader
                                        #   in Loop: Header=BB3_8 Depth=1
	li	t5, 0
	srliw	t4, t0, 1
	slli	t2, a6, 8
	add	t2, t2, a4
	slli	t3, t4, 2
.LBB3_11:                               # %.lr.ph.split
                                        #   Parent Loop BB3_8 Depth=1
                                        # =>  This Loop Header: Depth=2
                                        #       Child Loop BB3_12 Depth 3
	mulw	a0, t0, t5
	slli	a0, a0, 2
	add	a1, a7, a0
	add	a0, a1, t3
	mv	a2, t4
	mv	s0, t2
.LBB3_12:                               #   Parent Loop BB3_8 Depth=1
                                        #     Parent Loop BB3_11 Depth=2
                                        # =>    This Inner Loop Header: Depth=3
	vsetvli	s1, a2, e32, m4, ta, mu
	vle32.v	v8, (a0)
	vle32.v	v12, (s0)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a3, v8
	vwmul.vx	v24, v20, t6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a3
	vmnot.m	v0, v8
	vnmsub.vx	v16, a3, v20
	vmslt.vx	v8, v16, a3
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a5, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a5, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a0)
	sub	a2, a2, s1
	slli	s1, s1, 2
	add	a1, a1, s1
	add	a0, a0, s1
	add	s0, s0, s1
	bnez	a2, .LBB3_12
# %bb.13:                               #   in Loop: Header=BB3_11 Depth=2
	addi	t5, t5, 1
	bne	t5, t1, .LBB3_11
# %bb.14:                               # %._crit_edge
                                        #   in Loop: Header=BB3_8 Depth=1
	slli	t0, t0, 1
	addi	a0, a6, -1
	bnez	a6, .LBB3_8
# %bb.15:
	ld	ra, 24(sp)                      # 8-byte Folded Reload
	ld	s0, 16(sp)                      # 8-byte Folded Reload
	ld	s1, 8(sp)                       # 8-byte Folded Reload
	addi	sp, sp, 32
	ret
.LBB3_16:
	lui	a0, %hi(.L.str)
	addi	a0, a0, %lo(.L.str)
	lui	a1, %hi(.L__func__.rvv_ntt_transform_asm_helper)
	addi	a2, a1, %lo(.L__func__.rvv_ntt_transform_asm_helper)
	lui	a1, %hi(.L.str.5)
	addi	a3, a1, %lo(.L.str.5)
	li	a1, 584
	call	__assert_func
.LBB3_17:
	lui	a0, %hi(.L.str)
	addi	a0, a0, %lo(.L.str)
	lui	a1, %hi(.L__func__.rvv_ntt_transform_asm_helper)
	addi	a2, a1, %lo(.L__func__.rvv_ntt_transform_asm_helper)
	lui	a1, %hi(.L.str.12)
	addi	a3, a1, %lo(.L.str.12)
	li	a1, 667
	call	__assert_func
.LBB3_18:
	lui	a0, %hi(.L.str)
	addi	a0, a0, %lo(.L.str)
	lui	a1, %hi(.L__func__.rvv_ntt_transform_asm_helper)
	addi	a2, a1, %lo(.L__func__.rvv_ntt_transform_asm_helper)
	lui	a1, %hi(.L.str.13)
	addi	a3, a1, %lo(.L.str.13)
	li	a1, 668
	call	__assert_func
.Lfunc_end3:
	.size	rvv_ntt_transform_asm_helper, .Lfunc_end3-rvv_ntt_transform_asm_helper
                                        # -- End function

	.type	.L.str,@object                  # @.str
	.section	.rodata.str1.1,"aMS",@progbits,1
.L.str:
	.asciz	"poly_mult_rvv.c"
	.size	.L.str, 16

	.type	.L.str.3,@object                # @.str.3
.L.str.3:
	.asciz	"%d "
	.size	.L.str.3, 4

	.type	.L__func__.rvv_ntt_transform_asm_helper,@object # @__func__.rvv_ntt_transform_asm_helper
.L__func__.rvv_ntt_transform_asm_helper:
	.asciz	"rvv_ntt_transform_asm_helper"
	.size	.L__func__.rvv_ntt_transform_asm_helper, 30

	.type	.L.str.5,@object                # @.str.5
.L.str.5:
	.asciz	"_n > 1"
	.size	.L.str.5, 7

	.type	.L.str.7,@object                # @.str.7
.L.str.7:
	.asciz	"n == 2"
	.size	.L.str.7, 7

	.type	.L.str.8,@object                # @.str.8
.L.str.8:
	.asciz	"n == 2 && local_level == 6"
	.size	.L.str.8, 27

	.type	.L.str.9,@object                # @.str.9
.L.str.9:
	.asciz	"n == 4 && local_level == 5"
	.size	.L.str.9, 27

	.type	.L.str.10,@object               # @.str.10
.L.str.10:
	.asciz	"n <= FUNC_LMUL(__riscv_vsetvlmax_e32)()"
	.size	.L.str.10, 40

	.type	.L.str.11,@object               # @.str.11
.L.str.11:
	.asciz	"n == 8 && local_level == 4"
	.size	.L.str.11, 27

	.type	.L__func__.rvv_ntt_transform_fastest_helper,@object # @__func__.rvv_ntt_transform_fastest_helper
.L__func__.rvv_ntt_transform_fastest_helper:
	.asciz	"rvv_ntt_transform_fastest_helper"
	.size	.L__func__.rvv_ntt_transform_fastest_helper, 33

	.type	.L.str.12,@object               # @.str.12
.L.str.12:
	.asciz	"4 <= FUNC_LMUL(__riscv_vsetvlmax_e32)()"
	.size	.L.str.12, 40

	.type	.L.str.13,@object               # @.str.13
.L.str.13:
	.asciz	"8 <= FUNC_LMUL(__riscv_vsetvlmax_e32)()"
	.size	.L.str.13, 40

	.type	.L__func__.rvv_ntt_transform,@object # @__func__.rvv_ntt_transform
.L__func__.rvv_ntt_transform:
	.asciz	"rvv_ntt_transform"
	.size	.L__func__.rvv_ntt_transform, 18