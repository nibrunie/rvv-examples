#include <riscv_vector.h>
#include <stddef.h>

/** transpose of a n x n matrix
 *
 * @param dst address of destination matrix
 * @param src address of source matrix
 * @param n matrix dimensions
 */
void matrix_transpose(float *dst,
                      float *src,
                      size_t n) 
{
    size_t i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j) dst[i * n + j] = src[j * n + i];
};

/** 4x4 matrix transpose using strided stores */
void matrix_transpose_intrinsics_4x4(float *dst,
                                     float *src) 
{
    unsigned i;
    for (i = 0; i < 4; ++i) {
        vfloat32m1_t row = __riscv_vle32_v_f32m1(src + 4 * i, 4);
        __riscv_vsse32(dst + i, sizeof(float) * 4,row, 4);
    }
};


/** n x n matrix transpose using strided stores */
void matrix_transpose_intrinsics(float *dst,
                                 float *src,
                                 size_t n) 
{
    for (size_t row_id = 0; row_id < n; ++row_id) { // input row-index
        size_t avl = n;
        float* row_src = src + row_id * n;
        float* row_dst = dst + row_id;
        while (avl > 0) {
            size_t vl = __riscv_vsetvl_e32m1(avl);
            vfloat32m1_t row = __riscv_vle32_v_f32m1(row_src, vl);
            __riscv_vsse32(row_dst, sizeof(float) * n, row, vl);
            avl -= vl;
            row_src += vl;
            row_dst += vl * n;
        }
    }
};


/** n x n matrix transpose using strided loads */
void matrix_transpose_intrinsics_loads(float *dst,
                                       float *src,
                                       size_t n) 
{
    for (size_t row_id = 0; row_id < n; ++row_id) { // input row-index
        size_t avl = n;
        float* col_src = src + row_id;
        float* row_dst = dst + row_id * n;
        while (avl > 0) {
            size_t vl = __riscv_vsetvl_e32m1(avl);
            vfloat32m1_t row = __riscv_vlse32_v_f32m1(col_src, sizeof(float) * n, vl);
            __riscv_vse32(row_dst, row, vl);
            avl -= vl;
            col_src += vl * n;
            row_dst += vl;
        }
    }
};