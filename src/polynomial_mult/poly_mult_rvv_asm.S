// void rvv_ntt_transform_asm_helper(ntt_t* dst, int* coeffs, int _n, int level, int rootPowers[8][64]); 
// a0: destination buffer
// a1: input buffer
// a2: number of elements (should be 128)
// a3: level (should be )
// a4: address for root powers 2D array
.globl	rvv_ntt_transform_asm_helper # -- Begin function rvv_ntt_transform_asm_helper
	.p2align	1
	.type	rvv_ntt_transform_asm_helper,@function
rvv_ntt_transform_asm_helper:       # @rvv_ntt_transform_asm_helper
# %bb.0:
	vsetvli	a6, zero, e32, m1, ta, ma
	ld	t0, 8(a0)
	beqz	a2, .LBB3_4
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
# %bb.3:                                # %._crit_edge.loopexit
	ld	t0, 8(a0)
.LBB3_4:                                # %._crit_edge
	vsetvli	a0, zero, e32, m4, ta, ma
	beqz	a2, .LBB3_7
# %bb.5:
	li	a0, -86
	vsetvli	zero, a6, e8, m1, ta, ma
	vmv.v.x	v8, a0
	vmv.v.i	v9, -16
	addi	a0, a4, 1024
	vsetvli	a1, zero, e32, m4, ta, ma
	vle32.v	v12, (a0)
	li	a0, -52
	vsetvli	zero, a6, e8, m1, ta, ma
	vmv.v.x	v10, a0
	addi	a0, a4, 1280
	vsetvli	a1, zero, e32, m4, ta, ma
	vle32.v	v16, (a0)
	li	a6, 32
	li	a7, -1
	lui	a0, 1
	addi	t1, a0, 943
	addi	a0, a0, -767
	lui	a1, 1048575
	addi	t2, a1, 767
	mv	a3, t0
.LBB3_6:                                # %.lr.ph103
                                        # =>This Inner Loop Header: Depth=1
	vsetvli	a1, a2, e32, m4, ta, ma
	vle32.v	v20, (a3)
	srli	a5, a1, 1
	vsetvli	zero, a5, e32, m2, ta, ma
	vnsrl.wx	v24, v20, a6
	vnsrl.wi	v26, v20, 0
	vwaddu.vv	v28, v24, v26
	vwmaccu.vx	v28, a7, v26
	vsetvli	zero, a1, e32, m4, ta, mu
	vmv1r.v	v0, v8
	vrsub.vi	v20, v20, 0, v0.t
	vadd.vv	v20, v20, v28
	vmv1r.v	v0, v10
	vmul.vv	v20, v20, v16, v0.t
	vslidedown.vi	v24, v20, 2
	vslideup.vi	v24, v20, 2, v0.t
	vrsub.vi	v20, v20, 0, v0.t
	vadd.vv	v24, v20, v24
	vwmul.vx	v0, v24, t1
	vnsra.wi	v20, v0, 24
	vnmsub.vx	v20, a0, v24
	vmslt.vx	v11, v20, a0
	vmnot.m	v0, v11
	vadd.vx	v20, v20, t2, v0.t
	vmv1r.v	v0, v9
	vmul.vv	v20, v20, v12, v0.t
	vslidedown.vi	v24, v20, 4
	vslideup.vi	v24, v20, 4, v0.t
	vrsub.vi	v20, v20, 0, v0.t
	vadd.vv	v20, v20, v24
	vwmul.vx	v24, v20, t1
	vnsra.wi	v4, v24, 24
	vnmsub.vx	v4, a0, v20
	vmslt.vx	v11, v4, a0
	vmnot.m	v0, v11
	vadd.vx	v4, v4, t2, v0.t
	vse32.v	v4, (a3)
	sub	a2, a2, a1
	slli	a1, a1, 2
	add	a3, a3, a1
	bnez	a2, .LBB3_6
.LBB3_7:                                # %.lr.ph110
	addi	a6, a4, 768
	addi	a1, t0, 32
	li	a2, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	t1, a0, 767
	mv	a3, t0
	mv	a0, a6
.LBB3_8:                                # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a2, e32, m4, ta, mu
	vle32.v	v8, (a1)
	vle32.v	v12, (a0)
	vle32.v	v16, (a3)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a3)
	vse32.v	v12, (a1)
	sub	a2, a2, a5
	slli	a5, a5, 2
	add	a3, a3, a5
	add	a1, a1, a5
	add	a0, a0, a5
	bnez	a2, .LBB3_8
# %bb.9:
	addi	a1, t0, 64
	addi	a2, t0, 96
	li	a3, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a6
.LBB3_10:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	a5, a5, a0
	bnez	a3, .LBB3_10
# %bb.11:
	addi	a1, t0, 128
	addi	a2, t0, 160
	li	a3, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a6
.LBB3_12:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	a5, a5, a0
	bnez	a3, .LBB3_12
# %bb.13:
	addi	a1, t0, 192
	addi	a2, t0, 224
	li	a3, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a6
.LBB3_14:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	a5, a5, a0
	bnez	a3, .LBB3_14
# %bb.15:
	addi	a1, t0, 256
	addi	a2, t0, 288
	li	a3, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a6
.LBB3_16:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	a5, a5, a0
	bnez	a3, .LBB3_16
# %bb.17:
	addi	a1, t0, 320
	addi	a2, t0, 352
	li	a3, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a6
.LBB3_18:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	a5, a5, a0
	bnez	a3, .LBB3_18
# %bb.19:
	addi	a1, t0, 384
	addi	a2, t0, 416
	li	a3, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a6
.LBB3_20:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	a5, a5, a0
	bnez	a3, .LBB3_20
# %bb.21:
	addi	a1, t0, 448
	addi	a2, t0, 480
	li	a3, 8
	lui	a0, 1
	addi	a7, a0, 943
	addi	a0, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
.LBB3_22:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a6)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a0, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a0
	vmnot.m	v0, v8
	vnmsub.vx	v16, a0, v20
	vmslt.vx	v8, v16, a0
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t1, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t1, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a5
	slli	a5, a5, 2
	add	a1, a1, a5
	add	a2, a2, a5
	add	a6, a6, a5
	bnez	a3, .LBB3_22
# %bb.23:                               # %.lr.ph110.1
	addi	t1, a4, 512
	addi	a1, t0, 64
	li	a2, 16
	lui	a0, 1
	addi	a6, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	a7, a0, 767
	mv	a3, t0
	mv	a0, t1
.LBB3_24:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a2, e32, m4, ta, mu
	vle32.v	v8, (a1)
	vle32.v	v12, (a0)
	vle32.v	v16, (a3)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a7, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a7, v0.t
	vse32.v	v16, (a3)
	vse32.v	v12, (a1)
	sub	a2, a2, a5
	slli	a5, a5, 2
	add	a3, a3, a5
	add	a1, a1, a5
	add	a0, a0, a5
	bnez	a2, .LBB3_24
# %bb.25:
	addi	a1, t0, 128
	addi	a2, t0, 192
	li	a3, 16
	lui	a0, 1
	addi	a6, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	a7, a0, 767
	mv	a0, t1
.LBB3_26:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a0)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a7, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a7, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a5
	slli	a5, a5, 2
	add	a1, a1, a5
	add	a2, a2, a5
	add	a0, a0, a5
	bnez	a3, .LBB3_26
# %bb.27:
	addi	a1, t0, 256
	addi	a2, t0, 320
	li	a3, 16
	lui	a0, 1
	addi	a6, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	a7, a0, 767
	mv	a0, t1
.LBB3_28:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a0)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a7, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a7, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a5
	slli	a5, a5, 2
	add	a1, a1, a5
	add	a2, a2, a5
	add	a0, a0, a5
	bnez	a3, .LBB3_28
# %bb.29:
	addi	a1, t0, 384
	addi	a2, t0, 448
	li	a3, 16
	lui	a0, 1
	addi	a6, a0, 943
	addi	a5, a0, -767
	lui	a0, 1048575
	addi	a7, a0, 767
.LBB3_30:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (t1)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a5, v8
	vwmul.vx	v24, v20, a6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a5
	vmnot.m	v0, v8
	vnmsub.vx	v16, a5, v20
	vmslt.vx	v8, v16, a5
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a7, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a7, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	t1, t1, a0
	bnez	a3, .LBB3_30
# %bb.31:                               # %.lr.ph110.2
	addi	t1, a4, 256
	addi	a1, t0, 128
	li	a2, 32
	lui	a0, 1
	addi	a6, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	a7, a0, 767
	mv	a3, t0
	mv	a0, t1
.LBB3_32:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a2, e32, m4, ta, mu
	vle32.v	v8, (a1)
	vle32.v	v12, (a0)
	vle32.v	v16, (a3)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, a6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, t2
	vmnot.m	v0, v8
	vnmsub.vx	v16, t2, v20
	vmslt.vx	v8, v16, t2
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a7, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a7, v0.t
	vse32.v	v16, (a3)
	vse32.v	v12, (a1)
	sub	a2, a2, a5
	slli	a5, a5, 2
	add	a3, a3, a5
	add	a1, a1, a5
	add	a0, a0, a5
	bnez	a2, .LBB3_32
# %bb.33:
	addi	a1, t0, 256
	addi	a2, t0, 384
	li	a3, 32
	lui	a0, 1
	addi	a6, a0, 943
	addi	a5, a0, -767
	lui	a0, 1048575
	addi	a7, a0, 767
.LBB3_34:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (t1)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a5, v8
	vwmul.vx	v24, v20, a6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a5
	vmnot.m	v0, v8
	vnmsub.vx	v16, a5, v20
	vmslt.vx	v8, v16, a5
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a7, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a7, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	t1, t1, a0
	bnez	a3, .LBB3_34
# %bb.35:                               # %.lr.ph110.3
	addi	a0, t0, 256
	li	a1, 64
	lui	a3, 1
	addi	a6, a3, 943
	addi	a3, a3, -767
	lui	a5, 1048575
	addi	a5, a5, 767
.LBB3_36:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a2, a1, e32, m4, ta, mu
	vle32.v	v8, (a0)
	vle32.v	v12, (a4)
	vle32.v	v16, (t0)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a6
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a3, v8
	vwmul.vx	v24, v20, a6
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a3
	vmnot.m	v0, v8
	vnmsub.vx	v16, a3, v20
	vmslt.vx	v8, v16, a3
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a5, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a5, v0.t
	vse32.v	v16, (t0)
	vse32.v	v12, (a0)
	sub	a1, a1, a2
	slli	a2, a2, 2
	add	t0, t0, a2
	add	a0, a0, a2
	add	a4, a4, a2
	bnez	a1, .LBB3_36
# %bb.37:                               # %._crit_edge111.3
	ret
.Lfunc_end3:
	.size	rvv_ntt_transform_asm_helper, .Lfunc_end3-rvv_ntt_transform_asm_helper
                                        # -- End function