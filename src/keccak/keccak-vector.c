#include <stdint.h>
#include <riscv_vector.h>

#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

/*
Implementation by the Keccak Team, namely, Guido Bertoni, Joan Daemen,
Michaël Peeters, Gilles Van Assche and Ronny Van Keer,
hereby denoted as "the implementer".

For more information, feedback or questions, please refer to our website:
https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

/*
================================================================
The purpose of this source file is to demonstrate a readable and compact
implementation of all the Keccak instances approved in the FIPS 202 standard,
including the hash functions and the extendable-output functions (XOFs).

We focused on clarity and on source-code compactness,
rather than on the performance.

The advantages of this implementation are:
    + The source code is compact, after removing the comments, that is. :-)
    + There are no tables with arbitrary constants.
    + For clarity, the comments link the operations to the specifications using
        the same notation as much as possible.
    + There is no restriction in cryptographic features. In particular,
        the SHAKE128 and SHAKE256 XOFs can produce any output length.
    + The code does not use much RAM, as all operations are done in place.

The drawbacks of this implementation are:
    - There is no message queue. The whole message must be ready in a buffer.
    - It is not optimized for performance.

The implementation is even simpler on a little endian platform. Just define the
LITTLE_ENDIAN symbol in that case.

For a more complete set of implementations, please refer to
the Keccak Code Package at https://github.com/gvanas/KeccakCodePackage

For more information, please refer to:
    * [Keccak Reference] https://keccak.team/files/Keccak-reference-3.0.pdf
    * [Keccak Specifications Summary] https://keccak.team/keccak_specs_summary.html

This file uses UTF-8 encoding, as some comments use Greek letters.
================================================================
*/

/**
  * Function to compute the Keccak[r, c] sponge function over a given input.
  * @param  rate            The value of the rate r.
  * @param  capacity        The value of the capacity c.
  * @param  input           Pointer to the input message.
  * @param  inputByteLen    The number of input bytes provided in the input message.
  * @param  delimitedSuffix Bits that will be automatically appended to the end
  *                         of the input message, as in domain separation.
  *                         This is a byte containing from 0 to 7 bits
  *                         These <i>n</i> bits must be in the least significant bit positions
  *                         and must be delimited with a bit 1 at position <i>n</i>
  *                         (counting from 0=LSB to 7=MSB) and followed by bits 0
  *                         from position <i>n</i>+1 to position 7.
  *                         Some examples:
  *                             - If no bits are to be appended, then @a delimitedSuffix must be 0x01.
  *                             - If the 2-bit sequence 0,1 is to be appended (as for SHA3-*), @a delimitedSuffix must be 0x06.
  *                             - If the 4-bit sequence 1,1,1,1 is to be appended (as for SHAKE*), @a delimitedSuffix must be 0x1F.
  *                             - If the 7-bit sequence 1,1,0,1,0,0,0 is to be absorbed, @a delimitedSuffix must be 0x8B.
  * @param  output          Pointer to the buffer where to store the output.
  * @param  outputByteLen   The number of output bytes desired.
  * @pre    One must have r+c=1600 and the rate a multiple of 8 bits in this implementation.
  */
void Keccak(unsigned int rate, unsigned int capacity, const unsigned char *input, unsigned long long int inputByteLen, unsigned char delimitedSuffix, unsigned char *output, unsigned long long int outputByteLen);

/**
  *  Function to compute SHAKE128 on the input message with any output length.
  */
void FIPS202_SHAKE128(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen)
{
    Keccak(1344, 256, input, inputByteLen, 0x1F, output, outputByteLen);
}

/**
  *  Function to compute SHAKE256 on the input message with any output length.
  */
void FIPS202_SHAKE256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen)
{
    Keccak(1088, 512, input, inputByteLen, 0x1F, output, outputByteLen);
}

/**
  *  Function to compute SHA3-224 on the input message. The output length is fixed to 28 bytes.
  */
void FIPS202_SHA3_224(const unsigned char *input, unsigned int inputByteLen, unsigned char *output)
{
    Keccak(1152, 448, input, inputByteLen, 0x06, output, 28);
}

/**
  *  Function to compute SHA3-256 on the input message. The output length is fixed to 32 bytes.
  */
void FIPS202_SHA3_256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output)
{
    Keccak(1088, 512, input, inputByteLen, 0x06, output, 32);
}

/**
  *  Function to compute SHA3-384 on the input message. The output length is fixed to 48 bytes.
  */
void FIPS202_SHA3_384(const unsigned char *input, unsigned int inputByteLen, unsigned char *output)
{
    Keccak(832, 768, input, inputByteLen, 0x06, output, 48);
}

/**
  *  Function to compute SHA3-512 on the input message. The output length is fixed to 64 bytes.
  */
void FIPS202_SHA3_512(const unsigned char *input, unsigned int inputByteLen, unsigned char *output)
{
    Keccak(576, 1024, input, inputByteLen, 0x06, output, 64);
}

/*
================================================================
Technicalities
================================================================
*/

#include <stdint.h>

typedef uint64_t tKeccakLane;

#ifndef LITTLE_ENDIAN
/** Function to load a 64-bit value using the little-endian (LE) convention.
  * On a LE platform, this could be greatly simplified using a cast.
  */
static uint64_t load64(const uint8_t *x)
{
    int i;
    uint64_t u=0;

    for(i=7; i>=0; --i) {
        u <<= 8;
        u |= x[i];
    }
    return u;
}

/** Function to store a 64-bit value using the little-endian (LE) convention.
  * On a LE platform, this could be greatly simplified using a cast.
  */
static void store64(uint8_t *x, uint64_t u)
{
    unsigned int i;

    for(i=0; i<8; ++i) {
        x[i] = u;
        u >>= 8;
    }
}

/** Function to XOR into a 64-bit value using the little-endian (LE) convention.
  * On a LE platform, this could be greatly simplified using a cast.
  */
static void xor64(uint8_t *x, uint64_t u)
{
    unsigned int i;

    for(i=0; i<8; ++i) {
        x[i] ^= u;
        u >>= 8;
    }
}
#endif

/*
================================================================
A readable and compact implementation of the Keccak-f[1600] permutation.
================================================================
*/

#define ROL64(a, offset) ((((uint64_t)a) << offset) ^ (((uint64_t)a) >> (64-offset)))
#define i(x, y) ((x)+5*(y))

#ifdef LITTLE_ENDIAN
    #define readLane(x, y)          (((tKeccakLane*)state)[i(x, y)])
    #define writeLane(x, y, lane)   (((tKeccakLane*)state)[i(x, y)]) = (lane)
    #define XORLane(x, y, lane)     (((tKeccakLane*)state)[i(x, y)]) ^= (lane)
#else
    #define readLane(x, y)          load64((uint8_t*)state+sizeof(tKeccakLane)*i(x, y))
    #define writeLane(x, y, lane)   store64((uint8_t*)state+sizeof(tKeccakLane)*i(x, y), lane)
    #define XORLane(x, y, lane)     xor64((uint8_t*)state+sizeof(tKeccakLane)*i(x, y), lane)
#endif

/**
  * Function that computes the linear feedback shift register (LFSR) used to
  * define the round constants (see [Keccak Reference, Section 1.2]).
  */
int LFSR86540(uint8_t *LFSR)
{
    int result = ((*LFSR) & 0x01) != 0;
    if (((*LFSR) & 0x80) != 0)
        /* Primitive polynomial over GF(2): x^8+x^6+x^5+x^4+1 */
        (*LFSR) = ((*LFSR) << 1) ^ 0x71;
    else
        (*LFSR) <<= 1;
    return result;
}


/** This builtin was added recently to RVV intrinsics and seems to be missing from
 *  some recent clang version.
*/
vuint64m4_t __riscv_vcreate_v_u64m1_u64m4(vuint64m1_t v0, vuint64m1_t v1, vuint64m1_t v2, vuint64m1_t v3) {
    vuint64m4_t res = __riscv_vundefined_u64m4();
    res = __riscv_vset_v_u64m1_u64m4(res, 0, v0);
    res = __riscv_vset_v_u64m1_u64m4(res, 1, v1);
    res = __riscv_vset_v_u64m1_u64m4(res, 2, v2);
    res = __riscv_vset_v_u64m1_u64m4(res, 3, v3);
    return res;
}
/** This builtin was added recently to RVV intrinsics and seems to be missing from
 *  some recent clang version.
*/
vuint64m4_t __riscv_vrol_vx_u64m4(vuint64m4_t vs2, size_t rot, size_t vl){
    vuint64m4_t res = __riscv_vor_vv_u64m4(__riscv_vsll_vx_u64m4(vs2, rot,      vl),
                                           __riscv_vsrl_vx_u64m4(vs2, 64 - rot, vl),
                                           vl);
    return res;
}

vuint64m4_t __riscv_vrol_vv_u64m4(vuint64m4_t data, vuint64m4_t rots, size_t vl){
    vuint64m4_t rotsComp = __riscv_vrsub_vx_u64m4(rots, 64, vl); 
    vuint64m4_t res = __riscv_vor_vv_u64m4(__riscv_vsll_vv_u64m4(data, rots,      vl),
                                           __riscv_vsrl_vv_u64m4(data, rotsComp, vl),
                                           vl);
    return res;
}

/**
 * Function that computes the Keccak-f[1600] permutation on the given state.
 * original from: https://github.com/XKCP/XKCP/blob/master/Standalone/CompactFIPS202/C/Keccak-readable-and-compact.c
 */
void KeccakF1600_Round_vector(void *state, unsigned round, uint8_t* pLFSRstate)
{
    unsigned x, y, j, t;
    {   /* === θ step (see [Keccak Reference, Section 2.3.2]) === */
        tKeccakLane C[5], D;

        // state as a 5x5 matrix A[x, y] is stored with a column-major layout:
        // A[0,0]-A[1,0]-A[2,0]-A[3,0]-A[4,0] are contiguous in memory
        /* Compute the parity of the columns */
        /*
        tKeccakLane C_baseline[5];
        for(x=0; x<5; x++)
            C_baseline[x] = readLane(x, 0) ^ readLane(x, 1) ^ readLane(x, 2) ^ readLane(x, 3) ^ readLane(x, 4);
        //for(x=0; x<5; x++) {
        //    printf("C[%d]=%"PRIx64" vs C_rvv[%d]=%"PRIx64"\n", x, C[x], x, C_rvv[x]);
        //}
        // assert(!memcmp(C, C_baseline, 40));
        */
        /*
        vuint64m1x5_t rows_0 = __riscv_vlseg5e64_v_u64m1x5((uint64_t*)state, 2);
        vuint64m1_t row0_0 = __riscv_vget_v_u64m1x5_u64m1(rows_0, 0);
        vuint64m1_t row1_0 = __riscv_vget_v_u64m1x5_u64m1(rows_0, 1);
        vuint64m1_t row2_0 = __riscv_vget_v_u64m1x5_u64m1(rows_0, 2);
        vuint64m1_t row3_0 = __riscv_vget_v_u64m1x5_u64m1(rows_0, 3);
        vuint64m1_t row4_0 = __riscv_vget_v_u64m1x5_u64m1(rows_0, 4);

        vuint64m1x5_t rows_1 = __riscv_vlseg5e64_v_u64m1x5(((uint64_t*)state) + 10, 2);
        vuint64m1_t row0_1 = __riscv_vget_v_u64m1x5_u64m1(rows_1, 0);
        vuint64m1_t row1_1 = __riscv_vget_v_u64m1x5_u64m1(rows_1, 1);
        vuint64m1_t row2_1 = __riscv_vget_v_u64m1x5_u64m1(rows_1, 2);
        vuint64m1_t row3_1 = __riscv_vget_v_u64m1x5_u64m1(rows_1, 3);
        vuint64m1_t row4_1 = __riscv_vget_v_u64m1x5_u64m1(rows_1, 4);

        vuint64m1x5_t rows_2 = __riscv_vlseg5e64_v_u64m1x5(((uint64_t*)state) + 20, 1);
        vuint64m1_t row0_2 = __riscv_vget_v_u64m1x5_u64m1(rows_2, 0);
        vuint64m1_t row1_2 = __riscv_vget_v_u64m1x5_u64m1(rows_2, 1);
        vuint64m1_t row2_2 = __riscv_vget_v_u64m1x5_u64m1(rows_2, 2);
        vuint64m1_t row3_2 = __riscv_vget_v_u64m1x5_u64m1(rows_2, 3);
        vuint64m1_t row4_2 = __riscv_vget_v_u64m1x5_u64m1(rows_2, 4);

 
        vuint64m1_t zero_u64m1 = __riscv_vundefined_u64m1();
        vuint64m4_t row0 = __riscv_vcreate_v_u64m1_u64m4(row0_0, row0_1, row0_2, zero_u64m1);
        vuint64m4_t row1 = __riscv_vcreate_v_u64m1_u64m4(row1_0, row1_1, row1_2, zero_u64m1);
        vuint64m4_t row2 = __riscv_vcreate_v_u64m1_u64m4(row2_0, row2_1, row2_2, zero_u64m1);
        vuint64m4_t row3 = __riscv_vcreate_v_u64m1_u64m4(row3_0, row3_1, row3_2, zero_u64m1);
        vuint64m4_t row4 = __riscv_vcreate_v_u64m1_u64m4(row4_0, row4_1, row4_2, zero_u64m1);
        */
       vuint64m4_t row0 = __riscv_vle64_v_u64m4(((uint64_t*)state) + 0, 5);
       vuint64m4_t row1 = __riscv_vle64_v_u64m4(((uint64_t*)state) + 5, 5);
       vuint64m4_t row2 = __riscv_vle64_v_u64m4(((uint64_t*)state) + 10, 5);
       vuint64m4_t row3 = __riscv_vle64_v_u64m4(((uint64_t*)state) + 15, 5);
       vuint64m4_t row4 = __riscv_vle64_v_u64m4(((uint64_t*)state) + 20, 5);

        vuint64m4_t C_01 = __riscv_vxor_vv_u64m4(row0, row1, 5);
        vuint64m4_t C_23 = __riscv_vxor_vv_u64m4(row2, row3, 5);
        vuint64m4_t C_014 = __riscv_vxor_vv_u64m4(C_01, row4, 5);
        vuint64m4_t C_vector = __riscv_vxor_vv_u64m4(C_23, C_014, 5);

        __riscv_vse64_v_u64m4(C, C_vector, 5);
        /* Compute the θ effect for all the columns */
        uint64_t C_0 = __riscv_vmv_x_s_u64m4_u64(C_vector);
        vuint64m4_t C_4_ext = __riscv_vslidedown_vx_u64m4(C_vector, 4, 1);   // {C[4]}
        vuint64m4_t D_opLo = __riscv_vslide1down_vx_u64m4(C_vector, C_0, 5); // {C[1], C[2], C[3], C[4], C[0]}
        vuint64m4_t D_opHi = __riscv_vslide1up_vx_u64m4(C_vector, 0, 5);     // {0,    C[0], C[1], C[2], C[3]}
        D_opHi = __riscv_vxor_vv_u64m4_tu(D_opHi, D_opHi, C_4_ext, 1);

        // FIXME: intrinsics for vrol, currently not available in Docker's CLANG version 
        // D_opLo = __riscv_vrol_vx_u64m4(D_opLo, 1, 5);
        D_opLo = __riscv_vor_vv_u64m4(__riscv_vsll_vx_u64m4(D_opLo, 1,  5),
                                      __riscv_vsrl_vx_u64m4(D_opLo, 63, 5),
                                      5);
        vuint64m4_t D_rvv = __riscv_vxor_vv_u64m4(D_opLo, D_opHi, 5);
        tKeccakLane D_new[5];
        __riscv_vse64_v_u64m4(D_new, D_rvv, 5);
        
        /* Apply the θ effect */
        row0 = __riscv_vxor_vv_u64m4(row0, D_rvv, 5);
        row1 = __riscv_vxor_vv_u64m4(row1, D_rvv, 5);
        row2 = __riscv_vxor_vv_u64m4(row2, D_rvv, 5);
        row3 = __riscv_vxor_vv_u64m4(row3, D_rvv, 5);
        row4 = __riscv_vxor_vv_u64m4(row4, D_rvv, 5);

        __riscv_vse64_v_u64m4((uint64_t*)state,       row0, 5);
        __riscv_vse64_v_u64m4((uint64_t*)state + 5,  row1, 5);
        __riscv_vse64_v_u64m4((uint64_t*)state + 10,  row2, 5);
        __riscv_vse64_v_u64m4((uint64_t*)state + 15, row3, 5);
        __riscv_vse64_v_u64m4((uint64_t*)state + 20, row4, 5);

#if     0
        for(x=0; x<5; x++) {
            /* Compute the θ effect for a given column */
            D = C[(x+4)%5] ^ ROL64(C[(x+1)%5], 1);
            if (D != D_new[x]) {
                printf("D[%d] = %"PRIx64" vs D_new[%d]=%"PRIx64"\n", x, D, x, D_new[x]);
                assert(0);
            }
            /* Add the θ effect to the whole column */
            for (y=0; y<5; y++)
                XORLane(x, y, D);
        }
        assert(!memcmp(state, state_new, 200)); // , "state comparison after theta step");
#       endif
    }

    /* === ρ and π steps (see [Keccak Reference, Sections 2.3.3 and 2.3.4]) === */
    // indices and rotations generated with keccak/utils.py
    uint16_t offset_AtoB[] = {
        /*
        0, 6, 12, 18, 24, 
        3, 9, 10, 16, 22, 
        1, 7, 13, 19, 20, 
        4, 5, 11, 17, 23, 
        2, 8, 14, 15, 21, 
        */ 
        // byte offset for each index
        0, 48, 96, 144, 192, 
        24, 72, 80, 128, 176, 
        8, 56, 104, 152, 160, 
        32, 40, 88, 136, 184, 
        16, 64, 112, 120, 168, 
    };
    uint64_t rotation_B[] = {
        0, 44, 43, 21, 14, 
        28, 20, 3, 45, 61, 
        1, 6, 25, 8, 18, 
        27, 36, 10, 15, 56, 
        62, 55, 39, 41, 2, 
    };
    tKeccakLane B_rotated_rvv[25];

    unsigned rowId;
    // FIXME: could be done linearly on the 25-element array (no need to split by row)
    for (rowId = 0; rowId < 5; rowId++) {
        vuint16m1_t index_row = __riscv_vle16_v_u16m1(offset_AtoB + 5 * rowId, 5);
        vuint64m4_t B_row = __riscv_vluxei16_v_u64m4((uint64_t*)state, index_row, 5);
        vuint64m4_t rots_row = __riscv_vle64_v_u64m4(rotation_B + 5 * rowId, 5);
        B_row = __riscv_vrol_vv_u64m4(B_row, rots_row, 5);
        __riscv_vse64_v_u64m4(B_rotated_rvv + 5 * rowId, B_row, 5);
    }
    // using a second array to avoid overwriting elements
    // FIXME: to be optimized
    memcpy(state, B_rotated_rvv, 200);

#if 0
    {   /* === ρ and π steps (see [Keccak Reference, Sections 2.3.3 and 2.3.4]) === */
        tKeccakLane current, temp;
        /* Start at coordinates (1 0) */
        x = 1; y = 0;
        current = readLane(x, y);
        /* Iterate over ((0 1)(2 3))^t * (1 0) for 0 ≤ t ≤ 23 */
        for(t=0; t<24; t++) {
            /* Compute the rotation constant r = (t+1)(t+2)/2 */
            unsigned int r = ((t+1)*(t+2)/2)%64;
            /* Compute ((0 1)(2 3)) * (x y) */
            unsigned int Y = (2*x+3*y)%5; x = y; y = Y;
            /* Swap current and state(x,y), and rotate */
            temp = readLane(x, y);
            writeLane(x, y, ROL64(current, r));
            current = temp;
        }
    }
    assert(!memcmp(state, B_rotated_rvv, 200));
#endif

    {   /* === χ step (see [Keccak Reference, Section 2.3.1]) === */
        tKeccakLane temp[5];
        for(y=0; y<5; y++) {
            /* Take a copy of the plane */
            for(x=0; x<5; x++)
                temp[x] = readLane(x, y);
            /* Compute χ on the plane */
            for(x=0; x<5; x++)
                writeLane(x, y, temp[x] ^((~temp[(x+1)%5]) & temp[(x+2)%5]));
        }
    }

    {   /* === ι step (see [Keccak Reference, Section 2.3.5]) === */
        for(j=0; j<7; j++) {
            unsigned int bitPosition = (1<<j)-1; /* 2^j-1 */
            if (LFSR86540(pLFSRstate))
                XORLane(0, 0, (tKeccakLane)1<<bitPosition);
        }
    }

}

void KeccakF1600_StatePermute_vector(void *state)
{
    unsigned int round;
    uint8_t LFSRstate = 0x01;

    for(round=0; round<24; round++) {
        KeccakF1600_Round_vector(state, round, &LFSRstate);
    }
}

void KeccakF1600_StatePermute(void *state)
{
    KeccakF1600_StatePermute_vector(state);
}

/*
================================================================
A readable and compact implementation of the Keccak sponge functions
that use the Keccak-f[1600] permutation.
================================================================
*/

#define MIN(a, b) ((a) < (b) ? (a) : (b))

void Keccak(unsigned int rate, unsigned int capacity, const unsigned char *input, unsigned long long int inputByteLen, unsigned char delimitedSuffix, unsigned char *output, unsigned long long int outputByteLen)
{
    uint8_t state[200];
    unsigned int rateInBytes = rate/8;
    unsigned int blockSize = 0;
    unsigned int i;

    if (((rate + capacity) != 1600) || ((rate % 8) != 0))
        return;

    /* === Initialize the state === */
    memset(state, 0, sizeof(state));

    /* === Absorb all the input blocks === */
    while(inputByteLen > 0) {
        blockSize = MIN(inputByteLen, rateInBytes);
        for(i=0; i<blockSize; i++)
            state[i] ^= input[i];
        input += blockSize;
        inputByteLen -= blockSize;

        if (blockSize == rateInBytes) {
            KeccakF1600_StatePermute(state);
            blockSize = 0;
        }
    }

    /* === Do the padding and switch to the squeezing phase === */
    /* Absorb the last few bits and add the first bit of padding (which coincides with the delimiter in delimitedSuffix) */
    state[blockSize] ^= delimitedSuffix;
    /* If the first bit of padding is at position rate-1, we need a whole new block for the second bit of padding */
    if (((delimitedSuffix & 0x80) != 0) && (blockSize == (rateInBytes-1)))
        KeccakF1600_StatePermute(state);
    /* Add the second bit of padding */
    state[rateInBytes-1] ^= 0x80;
    /* Switch to the squeezing phase */
    KeccakF1600_StatePermute(state);

    /* === Squeeze out all the output blocks === */
    while(outputByteLen > 0) {
        blockSize = MIN(outputByteLen, rateInBytes);
        memcpy(output, state, blockSize);
        output += blockSize;
        outputByteLen -= blockSize;

        if (outputByteLen > 0)
            KeccakF1600_StatePermute(state);
    }
}