
	.globl	rvv_ntt_transform_asm_helper # -- Begin function rvv_ntt_transform_asm_helper
	.p2align	1
	.type	rvv_ntt_transform_asm_helper,@function
rvv_ntt_transform_asm_helper:       # @rvv_ntt_transform_asm_helper
# %bb.0:
	vsetvli	a7, zero, e32, m1, ta, ma
	ld	a6, 8(a0)
	beqz	a2, .LBB3_4
# %bb.1:                                # %.lr.ph.preheader
	lui	a3, %hi(ntt_coeff_indices_128)
	addi	a5, a3, %lo(ntt_coeff_indices_128)
	mv	a3, a2
.LBB3_2:                                # %.lr.ph
                                        # =>This Inner Loop Header: Depth=1
	vsetvli	t0, a3, e32, m4, ta, ma
	vle32.v	v8, (a5)
	vluxei32.v	v8, (a1), v8
	vse32.v	v8, (a6)
	sub	a3, a3, t0
	slli	t0, t0, 2
	add	a5, a5, t0
	add	a6, a6, t0
	bnez	a3, .LBB3_2
# %bb.3:                                # %._crit_edge.loopexit
	ld	a6, 8(a0)
.LBB3_4:                                # %._crit_edge
	vsetvli	a0, zero, e32, m4, ta, ma
	beqz	a2, .LBB3_7
# %bb.5:
	li	a0, -86
	vsetvli	zero, a7, e8, m1, ta, ma
	vmv.v.x	v8, a0
	vmv.v.i	v9, -16
	addi	a0, a4, 1024
	vsetvli	a1, zero, e32, m4, ta, ma
	vle32.v	v12, (a0)
	li	a0, -52
	vsetvli	zero, a7, e8, m1, ta, ma
	vmv.v.x	v10, a0
	addi	a0, a4, 1280
	vsetvli	a1, zero, e32, m4, ta, ma
	vle32.v	v16, (a0)
	lui	a1, 1
	addi	a7, a1, 943
	addi	a1, a1, -767
	lui	a3, 1048575
	addi	a5, a3, 767
	mv	a3, a6
.LBB3_6:                                # %.lr.ph101
                                        # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a2, e32, m4, ta, mu
	vle32.v	v20, (a3)
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
	vwmul.vx	v0, v24, a7
	vnsra.wi	v20, v0, 24
	vnmsub.vx	v20, a1, v24
	vmslt.vx	v11, v20, a1
	vmnot.m	v0, v11
	vadd.vx	v20, v20, a5, v0.t
	vmv1r.v	v0, v9
	vmul.vv	v20, v20, v12, v0.t
	vslidedown.vi	v24, v20, 4
	vslideup.vi	v24, v20, 4, v0.t
	vrsub.vi	v20, v20, 0, v0.t
	vadd.vv	v20, v20, v24
	vwmul.vx	v24, v20, a7
	vnsra.wi	v4, v24, 24
	vnmsub.vx	v4, a1, v20
	vmslt.vx	v11, v4, a1
	vmnot.m	v0, v11
	vadd.vx	v4, v4, a5, v0.t
	vse32.v	v4, (a3)
	sub	a2, a2, a0
	slli	a0, a0, 2
	add	a3, a3, a0
	bnez	a2, .LBB3_6
.LBB3_7:                                # %.lr.ph108
	addi	a7, a4, 768
	addi	a1, a6, 32
	li	a2, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	t1, a0, 767
	mv	a3, a6
	mv	a0, a7
.LBB3_8:                                # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a2, e32, m4, ta, mu
	vle32.v	v8, (a1)
	vle32.v	v12, (a0)
	vle32.v	v16, (a3)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, t0
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
	addi	a1, a6, 64
	addi	a2, a6, 96
	li	a3, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a7
.LBB3_10:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, t0
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
	addi	a1, a6, 128
	addi	a2, a6, 160
	li	a3, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a7
.LBB3_12:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, t0
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
	addi	a1, a6, 192
	addi	a2, a6, 224
	li	a3, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a7
.LBB3_14:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, t0
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
	addi	a1, a6, 256
	addi	a2, a6, 288
	li	a3, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a7
.LBB3_16:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, t0
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
	addi	a1, a6, 320
	addi	a2, a6, 352
	li	a3, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a7
.LBB3_18:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, t0
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
	addi	a1, a6, 384
	addi	a2, a6, 416
	li	a3, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	t2, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
	mv	a5, a7
.LBB3_20:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a5)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, t2, v8
	vwmul.vx	v24, v20, t0
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
	addi	a1, a6, 448
	addi	a2, a6, 480
	li	a3, 8
	lui	a0, 1
	addi	t0, a0, 943
	addi	a0, a0, -767
	lui	a5, 1048575
	addi	t1, a5, 767
.LBB3_22:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a7)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, t0
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a0, v8
	vwmul.vx	v24, v20, t0
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
	add	a7, a7, a5
	bnez	a3, .LBB3_22
# %bb.23:                               # %.lr.ph108.1
	addi	t1, a4, 512
	addi	a1, a6, 64
	li	a2, 16
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	t0, a0, 767
	mv	a3, a6
	mv	a0, t1
.LBB3_24:                               # =>This Inner Loop Header: Depth=1
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
	vadd.vx	v12, v12, t0, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t0, v0.t
	vse32.v	v16, (a3)
	vse32.v	v12, (a1)
	sub	a2, a2, a5
	slli	a5, a5, 2
	add	a3, a3, a5
	add	a1, a1, a5
	add	a0, a0, a5
	bnez	a2, .LBB3_24
# %bb.25:
	addi	a1, a6, 128
	addi	a2, a6, 192
	li	a3, 16
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	t0, a0, 767
	mv	a0, t1
.LBB3_26:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a0)
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
	vadd.vx	v12, v12, t0, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t0, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a5
	slli	a5, a5, 2
	add	a1, a1, a5
	add	a2, a2, a5
	add	a0, a0, a5
	bnez	a3, .LBB3_26
# %bb.27:
	addi	a1, a6, 256
	addi	a2, a6, 320
	li	a3, 16
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	t0, a0, 767
	mv	a0, t1
.LBB3_28:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a5, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (a0)
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
	vadd.vx	v12, v12, t0, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t0, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a5
	slli	a5, a5, 2
	add	a1, a1, a5
	add	a2, a2, a5
	add	a0, a0, a5
	bnez	a3, .LBB3_28
# %bb.29:
	addi	a1, a6, 384
	addi	a2, a6, 448
	li	a3, 16
	lui	a0, 1
	addi	a7, a0, 943
	addi	a5, a0, -767
	lui	a0, 1048575
	addi	t0, a0, 767
.LBB3_30:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (t1)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a5, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a5
	vmnot.m	v0, v8
	vnmsub.vx	v16, a5, v20
	vmslt.vx	v8, v16, a5
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t0, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t0, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	t1, t1, a0
	bnez	a3, .LBB3_30
# %bb.31:                               # %.lr.ph108.2
	addi	t1, a4, 256
	addi	a1, a6, 128
	li	a2, 32
	lui	a0, 1
	addi	a7, a0, 943
	addi	t2, a0, -767
	lui	a0, 1048575
	addi	t0, a0, 767
	mv	a3, a6
	mv	a0, t1
.LBB3_32:                               # =>This Inner Loop Header: Depth=1
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
	vadd.vx	v12, v12, t0, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t0, v0.t
	vse32.v	v16, (a3)
	vse32.v	v12, (a1)
	sub	a2, a2, a5
	slli	a5, a5, 2
	add	a3, a3, a5
	add	a1, a1, a5
	add	a0, a0, a5
	bnez	a2, .LBB3_32
# %bb.33:
	addi	a1, a6, 256
	addi	a2, a6, 384
	li	a3, 32
	lui	a0, 1
	addi	a7, a0, 943
	addi	a5, a0, -767
	lui	a0, 1048575
	addi	t0, a0, 767
.LBB3_34:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a0, a3, e32, m4, ta, mu
	vle32.v	v8, (a2)
	vle32.v	v12, (t1)
	vle32.v	v16, (a1)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a5, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a5
	vmnot.m	v0, v8
	vnmsub.vx	v16, a5, v20
	vmslt.vx	v8, v16, a5
	vmnot.m	v8, v8
	vadd.vx	v12, v12, t0, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, t0, v0.t
	vse32.v	v16, (a1)
	vse32.v	v12, (a2)
	sub	a3, a3, a0
	slli	a0, a0, 2
	add	a1, a1, a0
	add	a2, a2, a0
	add	t1, t1, a0
	bnez	a3, .LBB3_34
# %bb.35:                               # %.lr.ph108.3
	addi	a0, a6, 256
	li	a1, 64
	lui	a3, 1
	addi	a7, a3, 943
	addi	a3, a3, -767
	lui	a5, 1048575
	addi	a5, a5, 767
.LBB3_36:                               # =>This Inner Loop Header: Depth=1
	vsetvli	a2, a1, e32, m4, ta, mu
	vle32.v	v8, (a0)
	vle32.v	v12, (a4)
	vle32.v	v16, (a6)
	vmul.vv	v8, v8, v12
	vadd.vv	v20, v16, v8
	vsub.vv	v8, v16, v8
	vwmul.vx	v24, v8, a7
	vnsra.wi	v12, v24, 24
	vnmsub.vx	v12, a3, v8
	vwmul.vx	v24, v20, a7
	vnsra.wi	v16, v24, 24
	vmslt.vx	v8, v12, a3
	vmnot.m	v0, v8
	vnmsub.vx	v16, a3, v20
	vmslt.vx	v8, v16, a3
	vmnot.m	v8, v8
	vadd.vx	v12, v12, a5, v0.t
	vmv1r.v	v0, v8
	vadd.vx	v16, v16, a5, v0.t
	vse32.v	v16, (a6)
	vse32.v	v12, (a0)
	sub	a1, a1, a2
	slli	a2, a2, 2
	add	a6, a6, a2
	add	a0, a0, a2
	add	a4, a4, a2
	bnez	a1, .LBB3_36
# %bb.37:                               # %._crit_edge109.3
	ret
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