// Copyright 2022 Rivos Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "vlen-bits.h"

//
// Assembly routine signatures.
//

// VCLMUL

extern uint32_t
zvbc32e_vclmul_vv(
    uint32_t* dest,
    const uint32_t* src2,
    const uint32_t* src1,
    uint32_t n
);

extern uint32_t
zvbc32e_vclmul_vx(
    uint32_t* dest,
    const uint32_t* src2,
    uint32_t rs1,
    uint32_t n
);

// VCLMULH

extern uint32_t
zvbc32e_vclmulh_vv(
    uint32_t* dest,
    const uint32_t* src2,
    const uint32_t* src1,
    uint32_t n
);

extern uint32_t
zvbc32e_vclmulh_vx(
    uint32_t* dest,
    const uint32_t* src2,
    uint32_t rs1,
    uint32_t n
);

// SEW=16-bit VCLMUL

extern uint32_t
zvbc16e_vclmul_vv(
    uint16_t* dest,
    const uint16_t* src2,
    const uint16_t* src1,
    uint32_t n
);

extern uint32_t
zvbc16e_vclmul_vx(
    uint16_t* dest,
    const uint16_t* src2,
    uint16_t rs1,
    uint32_t n
);

// SEW=16-bit VCLMULH

extern uint32_t
zvbc16e_vclmulh_vv(
    uint16_t* dest,
    const uint16_t* src2,
    const uint16_t* src1,
    uint32_t n
);

extern uint32_t
zvbc16e_vclmulh_vx(
    uint16_t* dest,
    const uint16_t* src2,
    uint16_t rs1,
    uint32_t n
);

// SEW=8-bit VCLMUL

extern uint32_t
zvbc8e_vclmul_vv(
    uint8_t* dest,
    const uint8_t* src2,
    const uint8_t* src1,
    uint32_t n
);

extern uint32_t
zvbc8e_vclmul_vx(
    uint8_t* dest,
    const uint8_t* src2,
    uint8_t rs1,
    uint32_t n
);

// SEW=8-bit VCLMULH

extern uint32_t
zvbc8e_vclmulh_vv(
    uint8_t* dest,
    const uint8_t* src2,
    const uint8_t* src1,
    uint32_t n
);

extern uint32_t
zvbc8e_vclmulh_vx(
    uint8_t* dest,
    const uint8_t* src2,
    uint8_t rs1,
    uint32_t n
);


void
Assert(bool predicate, int lineno, const char* filename)
{
    if (predicate) {
        return;
    }
    fprintf(stderr, "\n%s: %d: Failed assertion.\n", filename, lineno);
    abort();
}

#define ASSERT(PREDICATE) Assert(PREDICATE, __LINE__, __FILE__)

typedef __uint128_t uint128_t;

// @brief Return a n-bit randomly generated number using rand()
//
// @return uint<n>_t
//
#define RAND_N(n) \
uint ##n##_t rand##n() { \
    return (uint##n##_t)rand(); \
}

RAND_N(32)
RAND_N(16)
RAND_N(8)

/** generic helper to compute a nxn carry-less multiply
 *  assuming 1 <= n <= 32.
 *  The full result (2n-1 bits) is returned.
 */
uint64_t
carryless_multiply(uint32_t a32, uint32_t b32, uint32_t n)
{
    const size_t kInputBits = n;
    const uint64_t b64 = b32;

    uint64_t accumulator = 0;
    for (size_t i = 0; i < kInputBits; ++i) {
        // If the i-th bit is set, then "add" (xor) a left-shifted
        // version of b to the accumulator.
        if (((a32 >> i) & 1) != 0) {
            accumulator ^= b64 << i;
        }
    }
    return accumulator;
}

uint64_t
carryless_multiply_32x32(uint32_t a32, uint32_t b32)
{
    return carryless_multiply(a32, b32, 32);
}



// @brief A helper function used to test the vclmul_vv instruction.
//
// This function inputs two vectors consisting of 64-bit scalars and outputs a
// vector such that c[i] is the high 64 bits of the 64x64->128 bit carryless
// multiplication of a[i] and b[i].
//
// @param a: The first vector we are multiplying
// @param b: The second vector we are multiplying
// @param n: Size of a, b and c vectors
// @param c: The output vector
//
void
carryless_multiply32vv(const uint32_t* a,
                            const uint32_t* b,
                            uint32_t shift,
                            size_t n,
                            uint32_t* c)
{
    for (size_t i = 0; i < n; ++i) {
        const uint64_t result64 = carryless_multiply_32x32(a[i], b[i]);
        c[i] = (uint32_t)(result64 >> shift);
    }
}

// @brief A helper function used to test the vclmul_vx/vclmul_vi instructions.
//
// This function inputs a vectors consisting of 64-bit scalars and a scalar.
// Outputs a vector such that c[i] is the high 64 bits of the 64x64->128 bit
// carryless multiplication of a[i] and b.
//
// @param a: The first vector we are multiplying
// @param b: The scalar we are multiplying
// @param n: Size of a, b and c vectors
// @param c: The output vector
//
void
carryless_multiply32vx(const uint32_t* a,
                            const uint32_t b,
                            uint32_t shift,
                            size_t n,
                            uint32_t* c)
{
    for (size_t i = 0; i < n; ++i) {
        const uint64_t result64 = carryless_multiply_32x32(a[i], b);
        c[i] = (uint32_t)(result64 >> shift);
    }
}

void
carryless_multiply16vv(const uint16_t* a,
                       const uint16_t* b,
                       uint16_t shift,
                       size_t n,
                       uint16_t* c)
{
    for (size_t i = 0; i < n; ++i) {
        const uint64_t result64 = carryless_multiply(a[i], b[i], 16);
        c[i] = (uint16_t)(result64 >> shift);
    }
}

void
carryless_multiply16vx(const uint16_t* a,
                       const uint16_t b,
                       uint16_t shift,
                       size_t n,
                       uint16_t* c)
{
    for (size_t i = 0; i < n; ++i) {
        const uint64_t result64 = carryless_multiply(a[i], b, 16);
        c[i] = (uint16_t)(result64 >> shift);
    }
}

void carryless_multiply8vv(const uint8_t* a,
                       const uint8_t* b,
                       uint8_t shift,
                       size_t n,
                       uint8_t* c)
{
    for (size_t i = 0; i < n; ++i) {
        const uint64_t result64 = carryless_multiply(a[i], b[i], 8);
        c[i] = (uint8_t)(result64 >> shift);
    }
}

void carryless_multiply8vx(const uint8_t* a,
                       const uint8_t b,
                       uint8_t shift,
                       size_t n,
                       uint8_t* c)
{
    for (size_t i = 0; i < n; ++i) {
        const uint64_t result64 = carryless_multiply(a[i], b, 8);
        c[i] = (uint8_t)(result64 >> shift);
    }
}

#define kNumElements 33
#define kRounds 1000

// @brief Tests vectorized carryless multiply intrinsic in assembly
// using randomly generated test vectors
//
// @return int 0 if intrinsics worked, 1 if they failed
//
int
test_rand_vclmul_32e(int high)
{

    uint32_t vs1[kNumElements];
    uint32_t vs2[kNumElements];

    uint32_t shift = high ? 32 : 0;

    if (high) {
        LOG("--- Testing Vectorized Carryless Multiply High (vector vector) 32-bit elements");
    } else {
        LOG("--- Testing Vectorized Carryless Multiply [Low] (vector vector) 32-bit elements");
    }

    for (size_t round = 0; round < kRounds; ++round) {
        for (size_t i = 0; i < kNumElements; ++i) {
            vs1[i] = rand32();
            vs2[i] = rand32();
        }
        uint32_t expected[kNumElements];
        carryless_multiply32vv(vs2, vs1, shift, kNumElements, expected);

        uint32_t actual[kNumElements] = { 0 };
        const uint32_t processed =
          high ? zvbc32e_vclmulh_vv(actual, vs2, vs1, kNumElements) : 
          zvbc32e_vclmul_vv(actual, vs2, vs1, kNumElements);

        if (processed != kNumElements ||
            memcmp(actual, expected, sizeof(actual))) {
            LOG("FAILURE: 'actual' does NOT match 'expected'");
            for (size_t i = 0; i < kNumElements; ++i) {
                LOG("expected[%3zd]: 0x%08" PRIx32
                    ", actual[%3zd]: 0x%08" PRIx32,
                    i, expected[i], i, actual[i]);
            }
            return 1;
        }
    }

    if (high) {
        LOG("--- Testing Vectorized Carryless Multiply High (vector scalar) 32-bit elements");
    } else {
        LOG("--- Testing Vectorized Carryless Multiply [Low] (vector scalar) 32-bit elements");
    }

    for (size_t round = 0; round < kRounds; ++round) {
        for (size_t i = 0; i < kNumElements; ++i) {
            vs2[i] = rand32();
        }
        const uint32_t rs1 = rand32();
        uint32_t expected[kNumElements];
        carryless_multiply32vx(vs2, rs1, shift, kNumElements, expected);

        uint32_t actual[kNumElements] = { 0 };
        const uint32_t processed = high ?
          zvbc32e_vclmulh_vx(actual, vs2, rs1, kNumElements) :
          zvbc32e_vclmul_vx(actual, vs2, rs1, kNumElements);

        if (processed != kNumElements ||
            memcmp(actual, expected, sizeof(actual))) {
            LOG("FAILURE: 'actual' does NOT match 'expected'");
            LOG(" - round    : %zu", round);
            LOG(" - processed: %" PRIu64, processed);
            for (size_t i = 0; i < kNumElements; ++i) {
                const uint32_t ac = actual[i];
                const uint32_t ex = expected[i];
                LOG("expected[%3zd]: 0x%08" PRIx32
                    ", actual[%3zd]: 0x%08" PRIx32 " %s",
                    i, ex, i, ac, (ac == ex ? "==" : "!="));
            }
            return 1;
        }
    }

    return 0;
}


// @brief Tests vectorized carryless multiply intrinsic in assembly
// using randomly generated test vectors
//
// @return int 0 if intrinsics worked, 1 if they failed
//
int
test_rand_vclmul_16e(int high)
{

    uint16_t vs1[kNumElements];
    uint16_t vs2[kNumElements];

    uint16_t shift = high ? 16 : 0;

    if (high) {
        LOG("--- Testing Vectorized Carryless Multiply High (vector vector) 16-bit elements");
    } else {
        LOG("--- Testing Vectorized Carryless Multiply [Low] (vector vector) 16-bit elements");
    }

    for (size_t round = 0; round < kRounds; ++round) {
        for (size_t i = 0; i < kNumElements; ++i) {
            vs1[i] = rand16();
            vs2[i] = rand16();
        }
        uint16_t expected[kNumElements];
        carryless_multiply16vv(vs2, vs1, shift, kNumElements, expected);

        uint16_t actual[kNumElements] = { 0 };
        const uint16_t processed =
          high ? zvbc16e_vclmulh_vv(actual, vs2, vs1, kNumElements) : 
          zvbc16e_vclmul_vv(actual, vs2, vs1, kNumElements);

        if (processed != kNumElements ||
            memcmp(actual, expected, sizeof(actual))) {
            LOG("FAILURE: 'actual' does NOT match 'expected'");
            for (size_t i = 0; i < kNumElements; ++i) {
                LOG("expected[%3zd]: 0x%04" PRIx16
                    ", actual[%3zd]: 0x%04" PRIx16,
                    i, expected[i], i, actual[i]);
            }
            return 1;
        }
    }

    if (high) {
        LOG("--- Testing Vectorized Carryless Multiply High (vector scalar) 16-bit elements");
    } else {
        LOG("--- Testing Vectorized Carryless Multiply [Low] (vector scalar) 16-bit elements");
    }

    for (size_t round = 0; round < kRounds; ++round) {
        for (size_t i = 0; i < kNumElements; ++i) {
            vs2[i] = rand16();
        }
        const uint16_t rs1 = rand16();
        uint16_t expected[kNumElements];
        carryless_multiply16vx(vs2, rs1, shift, kNumElements, expected);

        uint16_t actual[kNumElements] = { 0 };
        const uint16_t processed = high ?
          zvbc16e_vclmulh_vx(actual, vs2, rs1, kNumElements) :
          zvbc16e_vclmul_vx(actual, vs2, rs1, kNumElements);

        if (processed != kNumElements ||
            memcmp(actual, expected, sizeof(actual))) {
            LOG("FAILURE: 'actual' does NOT match 'expected'");
            LOG(" - round    : %zu", round);
            LOG(" - processed: %" PRIu64, processed);
            for (size_t i = 0; i < kNumElements; ++i) {
                const uint16_t ac = actual[i];
                const uint16_t ex = expected[i];
                LOG("expected[%3zd]: 0x%04" PRIx16
                    ", actual[%3zd]: 0x%04" PRIx16 " %s",
                    i, ex, i, ac, (ac == ex ? "==" : "!="));
            }
            return 1;
        }
    }

    return 0;
}


// @brief Tests vectorized carryless multiply intrinsic in assembly
// using randomly generated test vectors
//
// @return int 0 if intrinsics worked, 1 if they failed
//
int
test_rand_vclmul_8e(int high)
{

    uint8_t vs1[kNumElements];
    uint8_t vs2[kNumElements];

    uint8_t shift = high ? 8 : 0;

    if (high) {
        LOG("--- Testing Vectorized Carryless Multiply High (vector vector) 8-bit elements");
    } else {
        LOG("--- Testing Vectorized Carryless Multiply [Low] (vector vector) 8-bit elements");
    }

    for (size_t round = 0; round < kRounds; ++round) {
        for (size_t i = 0; i < kNumElements; ++i) {
            vs1[i] = rand8();
            vs2[i] = rand8();
        }
        uint8_t expected[kNumElements];
        carryless_multiply8vv(vs2, vs1, shift, kNumElements, expected);

        uint8_t actual[kNumElements] = { 0 };
        const uint8_t processed =
          high ? zvbc8e_vclmulh_vv(actual, vs2, vs1, kNumElements) : 
          zvbc8e_vclmul_vv(actual, vs2, vs1, kNumElements);

        if (processed != kNumElements ||
            memcmp(actual, expected, sizeof(actual))) {
            LOG("FAILURE: 'actual' does NOT match 'expected'");
            for (size_t i = 0; i < kNumElements; ++i) {
                LOG("expected[%3zd]: 0x%02" PRIx8
                    ", actual[%3zd]: 0x%02" PRIx8,
                    i, expected[i], i, actual[i]);
            }
            return 1;
        }
    }

    if (high) {
        LOG("--- Testing Vectorized Carryless Multiply High (vector scalar) 8-bit elements");
    } else {
        LOG("--- Testing Vectorized Carryless Multiply [Low] (vector scalar) 8-bit elements");
    }

    for (size_t round = 0; round < kRounds; ++round) {
        for (size_t i = 0; i < kNumElements; ++i) {
            vs2[i] = rand8();
        }
        const uint32_t rs1 = rand8();
        uint8_t expected[kNumElements];
        carryless_multiply8vx(vs2, rs1, shift, kNumElements, expected);

        uint8_t actual[kNumElements] = { 0 };
        const uint32_t processed = high ?
          zvbc8e_vclmulh_vx(actual, vs2, rs1, kNumElements) :
          zvbc8e_vclmul_vx(actual, vs2, rs1, kNumElements);

        if (processed != kNumElements ||
            memcmp(actual, expected, sizeof(actual))) {
            LOG("FAILURE: 'actual' does NOT match 'expected'");
            LOG(" - round    : %zu", round);
            LOG(" - processed: %" PRIu64, processed);
            for (size_t i = 0; i < kNumElements; ++i) {
                const uint32_t ac = actual[i];
                const uint32_t ex = expected[i];
                LOG("expected[%3zd]: 0x%02" PRIx8
                    ", actual[%3zd]: 0x%02" PRIx8 " %s",
                    i, ex, i, ac, (ac == ex ? "==" : "!="));
            }
            return 1;
        }
    }

    return 0;
}



// @brief Calls test functions for our intrinsics
//
// @return int
//
int
main()
{
    const uint64_t vlen = vlen_bits();
    LOG("VLEN = %" PRIu64, vlen);

    int res = 0;

    // vclmul SEW=32-bit
    res = test_rand_vclmul_32e(0 /* high=0 => low*/);
    if (res != 0) {
        return res;
    }

    // vclmulh SEW=32-bit
    res = test_rand_vclmul_32e(1 /* high */);
    if (res != 0) {
        return res;
    }

    // vclmul SEW=16-bit
    res = test_rand_vclmul_16e(0 /* high=0 => low*/);
    if (res != 0) {
        return res;
    }

    // vclmulh SEW=16-bit
    res = test_rand_vclmul_16e(1 /* high */);
    if (res != 0) {
        return res;
    }

    // vclmul SEW=8-bit
    res = test_rand_vclmul_8e(0 /* high=0 => low*/);
    if (res != 0) {
        return res;
    }

    // vclmulh SEW=8-bit
    res = test_rand_vclmul_8e(1 /* high */);
    if (res != 0) {
        return res;
    }

    return 0;
}
