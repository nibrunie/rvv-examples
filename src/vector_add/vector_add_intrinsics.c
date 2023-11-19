#include <riscv_vector.h>
#include <stddef.h>

void vector_add(float *dst,
                float *lhs,
                float *rhs,
                size_t avl)
{
    for (size_t vl; avl > 0; avl -= vl, lhs += vl, rhs += vl, dst += vl)
    {
        // compute loop body vector length from avl
        // (application vector length)
        vl = __riscv_vsetvl_e32m1(avl);
        // loading operands
        vfloat32m1_t vec_src_lhs = __riscv_vle32_v_f32m1(lhs, vl);
        vfloat32m1_t vec_src_rhs = __riscv_vle32_v_f32m1(rhs, vl);
        // actual vector addition
        vfloat32m1_t vec_acc = __riscv_vfadd_vv_f32m1(vec_src_lhs,
                                                      vec_src_rhs,
                                                      vl);
        // storing results
        __riscv_vse32_v_f32m1(dst, vec_acc, vl);
    }
}
