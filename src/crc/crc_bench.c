#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef MSG_SIZES
# define MSG_SIZES (1024 * 1024)
#endif 


uint32_t crc32_le_generic(uint32_t crc, unsigned char const *p,
					  size_t len,
					  uint32_t polynomial);

uint32_t crc32_be_generic(uint32_t crc, unsigned char const *p,
					  size_t len,
					  uint32_t polynomial);

void crc32init_le(const uint32_t polynomial);
void crc32init_be(const uint32_t polynomial);

/** RVV based implementation of Ethernet CRC32 */
uint32_t crcEth32_be_vector(uint32_t crc, unsigned char const *p, size_t len);
uint32_t crcEth32_le_vector(uint32_t crc, unsigned char const *p, size_t len);

/** RVV based optimized implementation of CRC32.
 *  It supports an arbitrary polynomial through the use of a table of reduction
 *  constants transmitted through a parameter.
 */
uint32_t crc_be_vector_opt(uint32_t crc, unsigned char const *p, size_t len, uint32_t redCsts[]);

/** basic counter read function */
unsigned long read_cycles(void)
{
  unsigned long cycles;
  asm volatile ("rdcycle %0" : "=r" (cycles));
  return cycles;
}

unsigned long read_instret(void)
{
  unsigned long instret;
  asm volatile ("rdinstret %0" : "=r" (instret));
  return instret;
}

int bench_crc_be(void)
{
    // Normal order, TODO: check if reverse is required
    uint32_t ethCRC32Poly = 0x04C11DB7;
    uint32_t start = 0, stop = 0;

    // initializing tables for CRC32 polynomial (baseline implementation)
    start = read_cycles();
    crc32init_be(ethCRC32Poly);
    stop = read_cycles();

#ifdef VERBOSE
    printf("CRC32 BE table init in %u cycle(s)\n", stop - start);
#endif // ifdef VERBOSE

    // generating constant for vector multiplication
    uint8_t dbgMsg[128] = {0};
    dbgMsg[0] = 1;
    // for (int i = 32; i <= 512; i +=32) {
    //    printf("CRC32(X^%d)  = %"PRIx32" (crc32_be_generic)\n", i, crc32_be_generic(0, dbgMsg, 1 + i / 8, ethCRC32Poly));
    //}

    start = read_cycles();
    uint32_t vectorRedConstants[2];
    vectorRedConstants[0] = crc32_be_generic(0, dbgMsg, 1 + 20, ethCRC32Poly); // X^160 mod polynomial
    vectorRedConstants[1] = crc32_be_generic(0, dbgMsg, 1 + 12, ethCRC32Poly); // X^160 mod polynomial
    stop = read_cycles();

#ifdef VERBOSE
    printf("crc_be_vector_opt reduction constant table init in %u cycle(s)\n", stop - start);
#endif // ifdef VERBOSE


    size_t msgSize = MSG_SIZES;
    unsigned char *inputMsg = (unsigned char*) malloc(msgSize);
    memset(inputMsg, 1, msgSize);

    start = read_cycles();
    // uint32_t resGeneric = crc32_le_generic(0, inputMsg, 1024, ethCRC32Poly);
    uint32_t resGeneric = crc32_be_generic(0, inputMsg, msgSize, ethCRC32Poly);
    stop = read_cycles();

    uint32_t delayGeneric = stop - start;
    float throughputGeneric = (double) delayGeneric / msgSize; 

    start = read_cycles();
    uint32_t resVector = crcEth32_be_vector(0, inputMsg, msgSize);
    stop = read_cycles();

    uint32_t delayVector = stop - start;
    float throughputVector = (double) delayVector / msgSize; 

    start = read_cycles();
    uint32_t resVectorOpt = crc_be_vector_opt(0, inputMsg, msgSize, vectorRedConstants);
    stop = read_cycles();

    uint32_t delayVectorOpt = stop - start;
    float throughputVectorOpt = (double) delayVectorOpt / msgSize; 

    // full message
#ifdef VERBOSE
    printf("CRC32(msg[%lu]) = %"PRIx32" (crc32_be_generic)   in %u cycle(s) [%.3f cycle(s) per Byte]\n", msgSize, resGeneric, delayGeneric, throughputGeneric);
    printf("CRC32(msg[%lu]) = %"PRIx32" (crcEth32_be_vector) in %u cycle(s) [%.3f cycle(s) per Byte]\n", msgSize, resVector,  delayVector, throughputVector);
    printf("CRC32(msg[%lu]) = %"PRIx32" (crc_be_vector_opt)  in %u cycle(s) [%.3f cycle(s) per Byte]\n", msgSize, resVectorOpt,  delayVectorOpt, throughputVectorOpt);
#else
    printf("BENCH %lu %"PRIx32" generic %u\n", msgSize, resGeneric, delayGeneric);
    printf("BENCH %lu %"PRIx32" _vector %u\n", msgSize, resVector,  delayVector);
    printf("BENCH %lu %"PRIx32" _vecopt %u\n", msgSize, resVectorOpt,  delayVectorOpt);
#endif // ifdef VERBOSE

    return !(resVector == resGeneric && resVectorOpt == resGeneric);
}

int bench_crc_le(void)
{
    // Normal order, TODO: check if reverse is required
    uint32_t ethCRC32PolyInv = 0xedb88320;
    uint32_t ethCRC32Poly    = 0x04C11DB7; // 0000 0100 1100 0001 0001 1101 1011 0111
    uint32_t start = 0, stop = 0;

    // initializing tables for CRC32 polynomial (baseline implementation)
    start = read_cycles();
    crc32init_le(ethCRC32PolyInv);
    stop = read_cycles();

#ifdef VERY_VERBOSE
    printf("CRC32 LE table init in %u cycle(s)\n", stop - start);
#endif

    // generating constant for vector multiplication
    uint8_t dbgMsg[128] = {0};
#ifdef VERY_VERBOSE
    dbgMsg[0] = 1;
    for (int i = 32; i <= 512; i +=32) {
        printf("CRC32(X^%d)  = %"PRIx32" (crc32_be_generic)\n", i, crc32_be_generic(0, dbgMsg, 1 + i / 8, ethCRC32Poly));
    }
    dbgMsg[0] = 0x80;
    for (int i = 32; i <= 128; i +=32) {
        printf("CRC32(X^%d)  = %"PRIx32" (crc32_le_generic)\n", i, crc32_le_generic(0, dbgMsg, 1 + i / 8, ethCRC32PolyInv));
    }
#endif

    start = read_cycles();
    uint32_t vectorRedConstants[2];
    vectorRedConstants[0] = crc32_le_generic(0, dbgMsg, 1 + 20, ethCRC32Poly); // X^160 mod polynomial
    vectorRedConstants[1] = crc32_le_generic(0, dbgMsg, 1 + 12, ethCRC32Poly); // X^160 mod polynomial
    stop = read_cycles();
    printf("crc_le_vector_opt reduction constant table init in %u cycle(s)\n", stop - start);


    size_t msgSize = MSG_SIZES;
    unsigned char *inputMsg = (unsigned char*) malloc(msgSize);
    memset(inputMsg, 1, msgSize);

    start = read_cycles();
    // uint32_t resGeneric = crc32_le_generic(0, inputMsg, 1024, ethCRC32Poly);
    uint32_t resGeneric = crc32_le_generic(0, inputMsg, msgSize, ethCRC32PolyInv);
    stop = read_cycles();

    uint32_t delayGeneric = stop - start;
    float throughputGeneric = (double) delayGeneric / msgSize; 

    start = read_cycles();
    uint32_t resVector = crcEth32_le_vector(0, inputMsg, msgSize);
    stop = read_cycles();

    uint32_t delayVector = stop - start;
    float throughputVector = (double) delayVector / msgSize; 

    // full message
#ifdef VERBOSE
    printf("CRC32(msg[%lu]) = %"PRIx32" (crc32_le_generic)   in %u cycle(s) [%.3f cycle(s) per Byte]\n", msgSize, resGeneric, delayGeneric, throughputGeneric);
    printf("CRC32(msg[%lu]) = %"PRIx32" (crcEth32_le_vector) in %u cycle(s) [%.3f cycle(s) per Byte]\n", msgSize, resVector,  delayVector, throughputVector);
#else
    printf("BENCH %lu %"PRIx32" generic %u\n", msgSize, resGeneric, delayGeneric);
    printf("BENCH %lu %"PRIx32" _vector %u\n", msgSize, resVector,  delayVector);
#endif // ifdef VERBOSE

    return !(resVector == resGeneric); //  && resVectorOpt == resGeneric);
}

int main(void) {
    return bench_crc_be() || bench_crc_le();
}