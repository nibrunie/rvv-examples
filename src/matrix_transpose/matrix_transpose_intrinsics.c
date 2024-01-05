#include <riscv_vector.h>
#include <stddef.h>

/** return the value of the instret counter
 *
 *  The instret counter counts the number of retired (executed) instructions.
*/
static unsigned long read_instret(void)
{
  unsigned long instret;
  asm volatile ("rdinstret %0" : "=r" (instret));
  return instret;
}

/** transpose of a n x n matrix
 *
 *  Baseline implementation (not using RVV explicitly through intrinsics)
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

/** Intrinsics based implementation of 4x4 32-bit matrix transpose */
void matrix_transpose_in_register(uint32_t* outputMat, uint32_t* inputMat) {
    vuint32m1_t matRegIn0 = __riscv_vle32_v_u32m1(inputMat, 4);
    vuint32m1_t matRegIn1 = __riscv_vle32_v_u32m1(inputMat + 4, 4);
    vuint32m1_t matRegIn2 = __riscv_vle32_v_u32m1(inputMat + 8, 4);
    vuint32m1_t matRegIn3 = __riscv_vle32_v_u32m1(inputMat + 12, 4);

    vbool32_t oddMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0xaaaa, 1));
    // vl=4 in the following
    // should be mapped to vslideup.vi
    vuint32m1_t smallTransposeMat0 = __riscv_vslideup_vx_u32m1_tumu(oddMask, matRegIn0, matRegIn1, 1, 4);
    vuint32m1_t smallTransposeMat2 = __riscv_vslideup_vx_u32m1_tumu(oddMask, matRegIn2, matRegIn3, 1, 4);

    vbool32_t evenMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0x5555, 1));
    // should me mapped to vslidedown.vi
    vuint32m1_t smallTransposeMat1 = __riscv_vslidedown_vx_u32m1_tumu(evenMask, matRegIn1, matRegIn0, 1, 4);
    vuint32m1_t smallTransposeMat3 = __riscv_vslidedown_vx_u32m1_tumu(evenMask, matRegIn3, matRegIn2, 1, 4);

    // should be mapped to vslideup.vi
    vuint32m1_t outMat0 = __riscv_vslideup_vx_u32m1_tu(smallTransposeMat0, smallTransposeMat2, 2, 4);
    vuint32m1_t outMat1 = __riscv_vslideup_vx_u32m1_tu(smallTransposeMat1, smallTransposeMat3, 2, 4);

    // vl=2 in the following
    // should me mapped to vslidedown.vi
    vuint32m1_t outMat2 = __riscv_vslidedown_vx_u32m1_tu(smallTransposeMat2, smallTransposeMat0, 2, 2);
    vuint32m1_t outMat3 = __riscv_vslidedown_vx_u32m1_tu(smallTransposeMat3, smallTransposeMat1, 2, 2);
    __riscv_vse32_v_u32m1(outputMat, outMat0, 4);
    __riscv_vse32_v_u32m1(outputMat + 4, outMat1, 4);
    __riscv_vse32_v_u32m1(outputMat + 8, outMat2, 4);
    __riscv_vse32_v_u32m1(outputMat + 12, outMat3, 4);
}


void matrix_transpose_segmented_load(float* dst, float* src) {
    asm volatile (
        // loading input matrix and performing the transpose at the same time
        // If the matrix was in registers, one could either use a unit-strided store
        // to store the input matrix in memory and then load it back with the 4-field
        // segmented load or one could use a 4-field segmented store to transpose
        // the matrix while storing it in a temporary memory buffer and then 
        // load it back (already transposed) in a vector register group through
        // the use of a unit-strided load
        "vsetivli a7, 4, e32, m1, tu, mu \n"
        "vlseg4e32.v v4, 0(%[inputMat])\n"
        "vsetivli a7, 16, e32, m4, tu, mu\n"
        "vse32.v v4, 0(%[outputMat])\n"
        :
        : [inputMat]"r"(src), [outputMat]"r"(dst)
        :
    );
}
void matrix_4x4_transpose_segmented_load_intrinsics(float* dst, float* src) {
    vfloat32m1x4_t data = __riscv_vlseg4e32_v_f32m1x4(src, 4);
    vfloat32m1_t data0 = __riscv_vget_v_f32m1x4_f32m1(data, 0);
    vfloat32m1_t data1 = __riscv_vget_v_f32m1x4_f32m1(data, 1);
    vfloat32m1_t data2 = __riscv_vget_v_f32m1x4_f32m1(data, 2);
    vfloat32m1_t data3 = __riscv_vget_v_f32m1x4_f32m1(data, 3);
    vfloat32m4_t packedData = __riscv_vcreate_v_f32m1_f32m4(data0, data1, data2, data3);
    __riscv_vse32_v_f32m4(dst, packedData, 16);
}

void matrix_4x4_transpose_segmented_store_intrinsics(float* dst, float* src) {
    vfloat32m4_t data = __riscv_vle32_v_f32m4(src, 16);
    vfloat32m1_t data0 = __riscv_vget_v_f32m4_f32m1(data, 0);
    vfloat32m1_t data1 = __riscv_vget_v_f32m4_f32m1(data, 1);
    vfloat32m1_t data2 = __riscv_vget_v_f32m4_f32m1(data, 2);
    vfloat32m1_t data3 = __riscv_vget_v_f32m4_f32m1(data, 3);
    vfloat32m1x4_t packedData = __riscv_vcreate_v_f32m1x4(data0, data1, data2, data3);
    __riscv_vsseg4e32_v_f32m1x4(dst, packedData, 4);
}

void matrix_transpose_segmented_store(float* dst, float* src) {
    asm volatile (
        "vsetivli a7, 16, e32, m4, tu, mu\n"
        "vle32.v v4, 0(%[inputMat])\n"
        "vsetivli a7, 4, e32, m1, tu, mu \n"
        "vsseg4e32.v v4, 0(%[outputMat])\n"
        :
        : [inputMat]"r"(src), [outputMat]"r"(dst)
        :
    );
}

/** extracting the array from the function matrix_4x4_transpose_vrgather
 *  seems to allow more optimization
*/
const uint16_t offsetIndex[] = {
    0, 4, 8, 12,
    1, 5, 9, 13,
    2, 6, 10, 14,
    3, 7, 11, 15
};

vfloat32m4_t matrix_4x4_transpose_vrgather(vfloat32m4_t src) {
    // loading array of indices
    vuint16m2_t indices = __riscv_vle16_v_u16m2(offsetIndex, 16);
    // performing transpose with in-register permutation
    vfloat32m4_t result = __riscv_vrgatherei16_vv_f32m4(src, indices, 16);
    return result;
}


/** To measure only the register-to-register count of instruction
 *  we wrap the bench to extract matrix load and store
*/
unsigned long matrix_4x4_transpose_vrgather_bench (float* dst, float* src) {
    unsigned long start, stop;
    vfloat32m4_t data = __riscv_vle32_v_f32m4(src, 16);
    start = read_instret();
    vfloat32m4_t result = matrix_4x4_transpose_vrgather(data);
    stop = read_instret();
    __riscv_vse32_v_f32m4(dst, result, 16);

    return stop - start;
}

vfloat32m4_t matrix_4x4_transpose_vslide(vfloat32m4_t src) {
    vfloat32m1_t matRegIn0 = __riscv_vget_v_f32m4_f32m1(src, 0);
    vfloat32m1_t matRegIn1 = __riscv_vget_v_f32m4_f32m1(src, 1);
    vfloat32m1_t matRegIn2 = __riscv_vget_v_f32m4_f32m1(src, 2);
    vfloat32m1_t matRegIn3 = __riscv_vget_v_f32m4_f32m1(src, 3);

    vbool32_t oddMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0xaaaa, 1));
    // vl=4 in the following
    // should be mapped to vslideup.vi
    vfloat32m1_t smallTransposeMat0 = __riscv_vslideup_vx_f32m1_tumu(oddMask, matRegIn0, matRegIn1, 1, 4);
    vfloat32m1_t smallTransposeMat2 = __riscv_vslideup_vx_f32m1_tumu(oddMask, matRegIn2, matRegIn3, 1, 4);

    vbool32_t evenMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0x5555, 1));
    // should me mapped to vslidedown.vi
    vfloat32m1_t smallTransposeMat1 = __riscv_vslidedown_vx_f32m1_tumu(evenMask, matRegIn1, matRegIn0, 1, 4);
    vfloat32m1_t smallTransposeMat3 = __riscv_vslidedown_vx_f32m1_tumu(evenMask, matRegIn3, matRegIn2, 1, 4);

    // should be mapped to vslideup.vi
    vfloat32m1_t outMat0 = __riscv_vslideup_vx_f32m1_tu(smallTransposeMat0, smallTransposeMat2, 2, 4);
    vfloat32m1_t outMat1 = __riscv_vslideup_vx_f32m1_tu(smallTransposeMat1, smallTransposeMat3, 2, 4);

    // vl=2 in the following
    // should me mapped to vslidedown.vi
    vfloat32m1_t outMat2 = __riscv_vslidedown_vx_f32m1_tu(smallTransposeMat2, smallTransposeMat0, 2, 2);
    vfloat32m1_t outMat3 = __riscv_vslidedown_vx_f32m1_tu(smallTransposeMat3, smallTransposeMat1, 2, 2);

    return __riscv_vcreate_v_f32m1_f32m4(outMat0, outMat1, outMat2, outMat3);
}

/** Intrinsics based implementation of 4x4 32-bit matrix transpose */
unsigned long matrix_4x4_transpose_vslide_bench(float* dst, float* src) {
    unsigned long start, stop;
    vfloat32m4_t data = __riscv_vle32_v_f32m4(src, 16);
    start = read_instret();
    vfloat32m4_t result = matrix_4x4_transpose_vslide(data);
    stop = read_instret();
    __riscv_vse32_v_f32m4(dst, result, 16);

    return stop - start;
}