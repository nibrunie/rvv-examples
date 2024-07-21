#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <riscv_vector.h>


typedef void (matrixMultiplyInterface_t)(uint32_t*, uint32_t*);

/** Intrinsics based implementation of 4x4 32-bit matrix transpose */
void transposeMatrix_4x4(uint32_t* outputMat, uint32_t* inputMat) {
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


/** Assembly version of RVV based 4x4 matrix transpose */
void transposeMatrix_4x4_asm(uint32_t* outputMat, uint32_t* inputMat) {
    // This routine implements a 4x4 transpose of a matrix of 32-bit elements.
    // It relies on vslideup/vslidedown, tail and masking.
    // It starts by doing 2x2 transpose between pairs of vector registers
    // and then builds the full 4x4 transpose by transposing tiles made from 2x2 blocks 
    //
    // The first step, 2x2 transpose is down through masked vslideup/vslidedown of 32-bit elements
    // The second step, 2x2 transpose of 2x2 tiles can be done through vslideup and either tail
    // undisturbed vslidedown of 32-bit or 64-bit elements or masked vslidedown.
    //
    // The code assumes vlen >= 128.
    asm volatile (
        // loading input matrix (this part should be removed if doing register-to-register transpose)
        "vsetivli a7, 16, e32, m4, tu, mu \n"
        "vle32.v v4, 0(%[inputMat])\n"

        "vsetivli a7, 4, e32, m1, tu, mu\n"
        // vbool32_t oddMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0xaaaa, 1));
        "vmv.v.i v0, 0xa\n" // materializing odd mask

        // vuint32m1_t smallTransposeMat0 = __riscv_vslideup_vx_u32m1_tumu(oddMask, matRegIn0, matRegIn1, 1, 4);
        "vmv.v.v v8, v4\n"       // v8 <- v4
        "vslideup.vi v8, v5, 1, v0.t\n"
        // vuint32m1_t smallTransposeMat2 = __riscv_vslideup_vx_u32m1_tumu(oddMask, matRegIn2, matRegIn3, 1, 4);
        "vmv.v.v v10, v6\n"       // v10 <- v6
        "vslideup.vi v10, v7, 1, v0.t\n"

        // vbool32_t evenMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0x5555, 1));
        "vmv.v.i v0, 0x5\n" // materializing even mask
        // vuint32m1_t smallTransposeMat1 = __riscv_vslidedown_vx_u32m1_tumu(evenMask, matRegIn1, matRegIn0, 1, 4);
        "vslidedown.vi v5, v4, 1, v0.t\n"
        // vuint32m1_t smallTransposeMat3 = __riscv_vslidedown_vx_u32m1_tumu(evenMask, matRegIn3, matRegIn2, 1, 4);
        "vslidedown.vi v7, v6, 1, v0.t\n"

        // vuint32m1_t outMat0 = __riscv_vslideup_vx_u32m1_tu(smallTransposeMat0, smallTransposeMat2, 2, 4);
        "vmv.v.v v4, v8\n"
        "vslideup.vi v4, v10, 2\n"
        // vuint32m1_t outMat1 = __riscv_vslideup_vx_u32m1_tu(smallTransposeMat1, smallTransposeMat3, 2, 4);
        "vslideup.vi v5, v7, 2\n"

        "vmv.v.v v6, v10\n"

        // to reduce code size the following fuseds two unit vslidedown into a single one with LMUL=2
        "vsetivli a7, 4, e64, m2, tu, mu\n"
        // vuint32m1_t outMat2 = __riscv_vslidedown_vx_u32m1_tu(smallTransposeMat2, smallTransposeMat0, 2, 2);
        // vuint32m1_t outMat3 = __riscv_vslidedown_vx_u32m1_tu(smallTransposeMat3, smallTransposeMat1, 2, 2);
        "vslidedown.vi v6, v8, 1, v0.t\n" // 2-vector group operation

        // storing result matrix
        "vsetivli a7, 16, e32, m4, tu, mu\n"
        "vse32.v v4, 0(%[outputMat])\n"
        :
        : [inputMat]"r"(inputMat), [outputMat]"r"(outputMat)
        :
    );
}

/** other assembly version of the 4x4 matrix transpose */
void transposeMatrix_4x4_asm_opt(uint32_t* outputMat, uint32_t* inputMat) {
    asm volatile (
        // loading input matrix
        "vsetivli a7, 16, e32, m4, tu, mu \n"
        "vle32.v v4, 0(%[inputMat])\n"

        "vsetivli a7, 16, e32, m4, tu, mu\n"
        "vmv.v.v v8, v4\n"       // v8v9v10v11 <- v4v5v6v7

        "vsetivli a7, 4, e32, m1, tu, mu\n"
        // vbool32_t oddMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0xaaaa, 1));
        "vmv.v.i v0, 0xa\n" // materializing odd mask

        // vuint32m1_t smallTransposeMat0 = __riscv_vslideup_vx_u32m1_tumu(oddMask, matRegIn0, matRegIn1, 1, 4);
        "vslideup.vi v8, v5, 1, v0.t\n"
        // vuint32m1_t smallTransposeMat2 = __riscv_vslideup_vx_u32m1_tumu(oddMask, matRegIn2, matRegIn3, 1, 4);
        "vslideup.vi v10, v7, 1, v0.t\n"

        // vbool32_t evenMask = __riscv_vreinterpret_v_u32m1_b32(__riscv_vmv_v_x_u32m1(0x5555, 1));
        "vmv.v.i v0, 0x5\n" // materializing even mask
        // vuint32m1_t smallTransposeMat1 = __riscv_vslidedown_vx_u32m1_tumu(evenMask, matRegIn1, matRegIn0, 1, 4);
        "vslidedown.vi v9, v4, 1, v0.t\n"
        // vuint32m1_t smallTransposeMat3 = __riscv_vslidedown_vx_u32m1_tumu(evenMask, matRegIn3, matRegIn2, 1, 4);
        "vslidedown.vi v11, v6, 1, v0.t\n"

        "vsetivli a7, 16, e32, m4, tu, mu\n"
        "vmv.v.v v4, v8\n"       // v8v9v10v11 -> v4v5v6v7

        "vsetivli a7, 4, e64, m2, tu, mu\n"
        // using already materialized mask 0x5
        // vuint32m1_t outMat2 = __riscv_vslidedown_vx_u32m1_tu(smallTransposeMat2, smallTransposeMat0, 2, 2);
        // vuint32m1_t outMat3 = __riscv_vslidedown_vx_u32m1_tu(smallTransposeMat3, smallTransposeMat1, 2, 2);
        "vslidedown.vi v6, v8, 1, v0.t\n"

        "vmv.v.i v0, 0xa\n" // materializing even mask
        // vuint32m1_t outMat0 = __riscv_vslideup_vx_u32m1_tu(smallTransposeMat0, smallTransposeMat2, 2, 4);
        // vuint32m1_t outMat1 = __riscv_vslideup_vx_u32m1_tu(smallTransposeMat1, smallTransposeMat3, 2, 4);
        "vslideup.vi v4, v10, 1, v0.t\n"


        // storing result matrix
        "vsetivli a7, 16, e32, m4, tu, mu\n"
        "vse32.v v4, 0(%[outputMat])\n"
        :
        : [inputMat]"r"(inputMat), [outputMat]"r"(outputMat)
        :
    );
}


void transposeMatrix_4x4_asm_mem(uint32_t* outputMat, uint32_t* inputMat) {
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
        : [inputMat]"r"(inputMat), [outputMat]"r"(outputMat)
        :
    );
}

/** basic counter read function */
unsigned long read_cycles(void)
{
  unsigned long cycles;
  asm volatile ("rdcycle %0" : "=r" (cycles));
  return cycles;
}

typedef struct {
    matrixMultiplyInterface_t *func;
    char label[10];
} matMulDesc_t;

int main(void) {
    // the input matrix is built such that its transpose should be a row-major matrix
    // of incremental coefficients {0, 1, 2, 3, 4, ...., 0xf}
    uint32_t inputMat[16] = {0, 4, 8, 0xc, 1, 5, 9, 0xd, 2, 6, 0xa, 0xe, 3, 7, 0xb, 0xf};
    uint32_t outputMat[16] = {0};

    matMulDesc_t descriptors[] = {
        {.func = transposeMatrix_4x4,         .label = "baseline"},
        {.func = transposeMatrix_4x4_asm,     .label = "asm"},
        {.func = transposeMatrix_4x4_asm_opt, .label = "asm_opt"},
        {.func = transposeMatrix_4x4_asm_mem, .label = "asm_mem"},
    };
    int i, j;
    printf("from:\n");
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) printf("| %d |", inputMat[i*4+j]);
        printf("\n");
    }

    int funcId;
    for (funcId = 0; funcId < sizeof(descriptors) / sizeof(matMulDesc_t); funcId++) {
        memset(outputMat, 64, 0);
        long start = read_cycles();
        descriptors[funcId].func(outputMat, inputMat);
        long stop = read_cycles();

        printf("to (%s) in %d cycles:\n", descriptors[funcId].label, stop - start);
        for (i = 0; i < 4; ++i) {
            for (j = 0; j < 4; ++j) printf("| %d |", outputMat[i*4+j]);
            printf("\n");
        }
    }

    // transposeMatrix_4x4(outputMat, inputMat);
    // transposeMatrix_4x4_asm(outputMatAsm, inputMat);
    return 0;
}