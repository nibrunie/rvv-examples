#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "bench_utils.h"

// setting a default message size (1MiB) if none is defined
#ifndef MSG_SIZES
# define MSG_SIZES (1024 * 1024)
#endif 


uint32_t crc32_le_generic(uint32_t crc, unsigned char const *p,
					  size_t len,
					  uint32_t polynomial);

uint32_t crc32_be_generic(uint32_t crc, unsigned char const *p,
					  size_t len,
					  uint32_t polynomial);

/** Initialization function (table init) for generic implementations */
void crc32init_le(const uint32_t polynomial);
void crc32init_be(const uint32_t polynomial);

const uint32_t ethCRC32Poly = 0x04C11DB7;
const uint32_t ethCRC32PolyInv = 0xedb88320;

// wrapper to build ETH32 specific versions of LE and BE generic CRC routines
uint32_t crcEth32_be_generic(uint32_t crc, unsigned char const *p, size_t len) {
    return crc32_be_generic(crc, p, len, ethCRC32Poly);
}
uint32_t crcEth32_le_generic(uint32_t crc, unsigned char const *p, size_t len) {
    return crc32_le_generic(crc, p, len, ethCRC32PolyInv);
}

/** RVV based implementation of Ethernet CRC32 */
uint32_t crcEth32_be_vector(uint32_t crc, unsigned char const *p, size_t len);
uint32_t crcEth32_le_vector(uint32_t crc, unsigned char const *p, size_t len);

/** RVV based optimized implementation of CRC32.
 *  It supports an arbitrary polynomial through the use of a table of reduction
 *  constants transmitted through a parameter.
 */
uint32_t crcEth32_be_vector_opt(uint32_t crc, unsigned char const *p, size_t len);

uint32_t crcEth32_le_vector_opt(uint32_t crc, unsigned char const *p, size_t len);

uint32_t crcEth32_be_vector_zvbc32e(uint32_t crc, unsigned char const *p, size_t len);


typedef struct {
    uint32_t (*crc_func)(uint32_t seed, unsigned char const* p, size_t len);
    bool be;
    char* label;
} crc_bench_t;

int bench_crc(void) {
    uint32_t start = 0, stop = 0;

    // initializing tables for CRC32 polynomial (baseline implementation)
    start = read_perf_counter();
    crc32init_le(ethCRC32PolyInv);
    stop = read_perf_counter();
#ifdef VERY_VERBOSE
    printf("CRC32 LE table init in %u " PERF_METRIC "(s)\n", stop - start);
#endif
    start = read_perf_counter();
    crc32init_be(ethCRC32Poly);
    stop = read_perf_counter();
#ifdef VERY_VERBOSE
    printf("CRC32 BE table init in %u " PERF_METRIC "(s)\n", stop - start);
#endif

    crc_bench_t benchmarks[] = {
        // labels should be padded to all have the same width (result display alignment)
        // CRC BE variants
        {.crc_func = crcEth32_be_generic,        .be = true, .label = "crcEth32_be_generic       "},
        {.crc_func = crcEth32_be_vector,         .be = true, .label = "crcEth32_be_vector        "},
        {.crc_func = crcEth32_be_vector_opt,     .be = true, .label = "crcEth32_be_vector_opt    "},
        {.crc_func = crcEth32_be_vector_zvbc32e, .be = true, .label = "crcEth32_be_vector_zvbc32e"},

        // CRC LE variants
        {.crc_func = crcEth32_le_generic,        .be = false, .label = "crcEth32_le_generic       "},
        {.crc_func = crcEth32_le_vector,         .be = false, .label = "crcEth32_le_vector        "},
        // FIXME: not developped yet
        // {.crc_func = crcEth32_le_vector_opt,     .be = false, .label = "crcEth32_le_vector_opt    "},
    };

    size_t msgSizes[] = {16, 32, 64, 128, 512, 1024, 16384, MSG_SIZES};
    int error = 0;
#   ifndef VERBOSE
        printf("message_size, crc_result, label, perf(" PERF_METRIC ")\n");
#   endif

    for (int size_id = 0; size_id < sizeof(msgSizes) / sizeof(size_t); ++size_id) {
        size_t msgSize = msgSizes[size_id];
        unsigned char *inputMsg = (unsigned char*) malloc(msgSize);
        memset(inputMsg, 1, msgSize);

        // computing golden
        uint32_t golden_be = crcEth32_be_generic(0, inputMsg, msgSize);
        uint32_t golden_le = crcEth32_le_generic(0, inputMsg, msgSize);


        for (int bench_id = 0; bench_id < sizeof(benchmarks) / sizeof(crc_bench_t); bench_id++) {
            start = read_perf_counter();
            uint32_t result = benchmarks[bench_id].crc_func(0, inputMsg, msgSize);
            stop = read_perf_counter();

            // checks
            error += result != (benchmarks[bench_id].be ? golden_be : golden_le);

            uint32_t delay = stop - start;
            float throughput = (double) delay / msgSize; 

            // full message
    #       ifdef VERBOSE
            printf("CRC32(msg[%lu]) = %"PRIx32" (%s)      in %u " PERF_METRIC "(s) [%.3f " PERF_METRIC "(s) per Byte]\n",
                msgSize, result, benchmarks[bench_id].label, delay, throughput);
    #       else
            printf("%lu, %"PRIx32", %s, %u\n", msgSize, result, benchmarks[bench_id].label, delay);
    #       endif

        }

        free(inputMsg);

    }

    return error;
}


void compute_crc_reduction_constants() {
    uint32_t start = 0, stop = 0;

    // generating constant for vector multiplication
    uint8_t dbgMsg[128] = {0};
    // LE CRC requires X^MSB to appear as the most significant bit of the first byte
    dbgMsg[0] = 0x80;
    for (int i = 32; i <= 128; i +=32) {
        printf("CRC32(X^%d)  = %"PRIx32" (crc32_le_generic)\n", i, crc32_le_generic(0, dbgMsg, 1 + i / 8, ethCRC32PolyInv));
    }
    start = read_perf_counter();
    uint32_t vectorRedConstantsLE[2];
    vectorRedConstantsLE[0] = crc32_le_generic(0, dbgMsg, 1 + 20, ethCRC32Poly); // X^160 mod polynomial
    vectorRedConstantsLE[1] = crc32_le_generic(0, dbgMsg, 1 + 12, ethCRC32Poly); // X^160 mod polynomial
    stop = read_perf_counter();
    printf("crc_le_vector_opt reduction constant table init in %u " PERF_METRIC "(s)\n", stop - start);
    for (int i = 0; i < 2; ++i) printf("[LE] vectorRedConstants[%d] = %"PRIx32"\n", i, vectorRedConstantsLE[i]);


    // BE CRC requires X^MSB to appear as the least significant bit of the first byte
    dbgMsg[0] = 1;
    for (int i = 32; i <= 512; i +=32) {
        printf("CRC32(X^%d)  = %"PRIx32" (crc32_be_generic)\n", i, crc32_be_generic(0, dbgMsg, 1 + i / 8, ethCRC32Poly));
    }

    start = read_perf_counter();
    uint32_t vectorRedConstantsBE[2];
    vectorRedConstantsBE[0] = crc32_be_generic(0, dbgMsg, 1 + 20, ethCRC32Poly); // X^160 mod polynomial
    vectorRedConstantsBE[1] = crc32_be_generic(0, dbgMsg, 1 + 12, ethCRC32Poly); // X^96 mod polynomial
    for (int i = 0; i < 2; ++i) printf("[BE] vectorRedConstants[%d] = %"PRIx32"\n", i, vectorRedConstantsBE[i]);

    uint32_t vectorRedConstantsZvbc32e[4];
    for (int i = 0; i < 4; ++i) {
        vectorRedConstantsZvbc32e[i] = crc32_be_generic(0, dbgMsg, 1 + 16 - 4 * i, ethCRC32Poly);
        printf("[BE] vectorRedConstantsZvbc32e[%d] = %"PRIx32"\n", i, vectorRedConstantsZvbc32e[i]);
    }
    stop = read_perf_counter();

#ifdef VERBOSE
    printf("crc_be_vector_opt reduction constant table init in %u " PERF_METRIC "(s)\n", stop - start);
#endif // ifdef VERBOSE
}


int main(void) {
    // define to compute and display the static constants (X^i mod P) = CRC(X^(i-32)) used
    // in folding implementations
#ifdef DISPLAY_REDUCTION_CSTS
    compute_crc_reduction_constants();
#endif
    return bench_crc();
}