#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <riscv_vector.h>

/** declaration of generic crc32_be_generic function, used to complete
 *  the vector implementations (epilog).
*/
uint32_t crc32_be_generic(uint32_t crc, unsigned char const *p,
					  size_t len,
					  uint32_t polynomial);

/** crc32Eth32_be_vector() - Calculate bitwise little-endian Ethernet AUTODIN II
 *			CRC32/CRC32C
 * @crc: seed value for computation.  ~0 for Ethernet, sometimes 0 for other
 *	 uses, or the previous crc32/crc32c value if computing incrementally.
 * @p: pointer to buffer over which CRC32/CRC32C is run
 * @len: length of buffer @p
 *
 */
uint32_t crcEth32_be_vector(uint32_t crc, unsigned char const *p, size_t len)
{
  uint32_t polynomial = 0x04C11DB7;
	int i;
  size_t avl = len / 4; // 4-byte per element 

  // we take a pair of 32-bit elements (E, F) from the message and reduce
  // E.X^32 by multiplying E.R where R=X^32[CRCPoly]
  // F.X^64 by multiplying F.S where S=X^64[CRCPoly]

   // {(X^96 mod P), (X^64 mod P)} = {0xe8a45605, 0xf200aa66};
   // {(X^64 mod P), (X^32 mod P)}
   const uint32_t redConstants[] = {0xf200aa66, 0x490d678d};
   vuint32mf2_t redConstantVector = __riscv_vle32_v_u32mf2(redConstants, 2);
   vuint64m1_t extRedCstVector = __riscv_vzext_vf2_u64m1(redConstantVector, 2);

   vuint64m1_t crcAcc = __riscv_vmv_v_x_u64m1(0, 2);
   const vuint64m1_t zeroVecU64M1 = __riscv_vmv_v_x_u64m1(0, 2);

  for (size_t vl; avl > 3; avl -= vl, p += 4 * vl, len -= 4*vl) {
      // compute loop body vector length from application vector length avl
      vl = __riscv_vsetvl_e32mf2(avl);
      // printf("vector loop body, vl=%lu, avl=%lu, len=%lu\n", vl, avl, len);
      // printf("p=%llx\n", p);
      vuint32mf2_t inputData = __riscv_vle32_v_u32mf2((uint32_t*) p, vl);
      // byte swapping each 32-bit element
      // vuint32m1_t __riscv_vreinterpret_v_u64m1_u32m1 (vuint64m1_t src);
      crcAcc = __riscv_vrev8_v_u64m1(crcAcc, 1); 
      // printf("post rev8 crcAccU64=%"PRIx64"\n", __riscv_vmv_x_s_u64m1_u64(crcAcc));
      vuint32m1_t crcAccU32 = __riscv_vreinterpret_v_u64m1_u32m1 (crcAcc);
      // printf("post rev8 crcAccU32=%"PRIx32"\n", __riscv_vmv_x_s_u32m1_u32(crcAccU32));
      // printf("post rev8 crcAccU32=%"PRIx32"\n", __riscv_vmv_x_s_u32m1_u32(__riscv_vslidedown_vx_u32m1(crcAccU32, 1, 2)));
      vuint32mf2_t crcAccU32mf2 = __riscv_vlmul_trunc_v_u32m1_u32mf2(crcAccU32);
      // printf("post rev8 crcAccU32mf2=%"PRIx32"\n", __riscv_vmv_x_s_u32mf2_u32(crcAccU32mf2));
      // printf("post rev8 crcAccU32mf2=%"PRIx32"\n", __riscv_vmv_x_s_u32mf2_u32(__riscv_vslidedown_vx_u32mf2(crcAccU32mf2, 1, 2)));
      // vuint32mf2_t __riscv_vlmul_trunc_v_u32m1_u32mf2 (vuint32m1_t op1);
      inputData = __riscv_vxor_vv_u32mf2(inputData, crcAccU32mf2, vl);
      inputData = __riscv_vrev8_v_u32mf2(inputData, vl);

      // expanding data u32 -> u64
#ifdef HAS_ZVBB_SUPPORT
      vuint64m1_t extInputData = __riscv_vwsll_vx_u64m1(inputData, 0, vl);
#else // No HAS_ZVBB_SUPPORT
      vuint64m1_t extInputData = __riscv_vzext_vf2_u64m1(inputData, vl);
#endif // HAS_ZVBB_SUPPORT
      //
      vuint64m1_t multRes = __riscv_vclmul_vv_u64m1(extInputData, extRedCstVector, vl);
      crcAcc = __riscv_vredxor_vs_u64m1_u64m1(multRes, zeroVecU64M1, vl);
      // printf("crcAccU64=%"PRIx64"\n", __riscv_vmv_x_s_u64m1_u64(crcAcc));
  }
  uint64_t crcAccBuffer[1] = {0};
  uint32_t preCrc = __riscv_vmv_x_s_u64m1_u64(crcAcc);
  // printf("preCrc=%"PRIx32"\n", preCrc);
  // printf("p=%llx\n", p);
  // byte order reversing to ensure MS-byte is stored first (lowest address)
  crcAcc = __riscv_vrev8_v_u64m1(crcAcc, 1);
  __riscv_vse64_v_u64m1(crcAccBuffer, crcAcc, 1);
  uint8_t* crcAccBufferU8 = (uint8_t*) crcAccBuffer;

  // printf("vector loop end len=%lu\n", len);
  const uint32_t ethCRC32Poly = 0x04C11DB7;
  return crc32_be_generic(0, crcAccBufferU8, 8, ethCRC32Poly) ^ crc32_be_generic(0, p, len, ethCRC32Poly);
}


/** crc32Eth32_be_vector() - Calculate bitwise little-endian Ethernet AUTODIN II
 *			CRC32/CRC32C
 * @crc: seed value for computation.  ~0 for Ethernet, sometimes 0 for other
 *	 uses, or the previous crc32/crc32c value if computing incrementally.
 * @p: pointer to buffer over which CRC32/CRC32C is run
 * @len: length of buffer @p
 * @redCsts: reduction consts {X^160, X^96, ...}
 *
 */
uint32_t crcEth32_be_vector_opt(uint32_t crc, unsigned char const *p, size_t len)
{
	int i;
  size_t avl = len / 8; // 8-byte per element

  // we take a pair of 64-bit elements (E, F) from the message and reduce
  // E.X^96 by multiplying  E.R where R=X^96[CRCPoly]
  // F.X^160 by multiplying F.S where S=X^160[CRCPoly]
  // {(X^160 mod P), (X^96 P)} are stored in redCsts table
  const uint32_t redCsts[2] = {
    0xc5b9cd4c,
    0xe8a45605,
  };

   vuint32mf2_t redConstantVector = __riscv_vle32_v_u32mf2(redCsts, 2);
   vuint64m1_t extRedCstVector = __riscv_vzext_vf2_u64m1(redConstantVector, 2);

   vuint64m1_t crcAccLo = __riscv_vmv_v_x_u64m1(0, 2);
   vuint64m1_t crcAccHi = __riscv_vmv_v_x_u64m1(0, 2);
   vuint64m1_t crcAcc = __riscv_vmv_v_x_u64m1(0, 2);
   const vuint64m1_t zeroVecU64M1 = __riscv_vmv_v_x_u64m1(0, 2);

   // hoisting setting vl (assuming it won't change in loop)
   size_t vl = __riscv_vsetvl_e64m1(-1);

  for (; avl >= 2*vl; avl -= vl, p += 8 * vl, len -= 8*vl) {
      // compute loop body vector length from application vector length avl
      assert(vl == 2);
      vuint64m1_t inputData = __riscv_vle64_v_u64m1((uint64_t*) p, vl);
      // byte swapping the data to align their endianess with the CRC accumulator
      inputData = __riscv_vrev8_v_u64m1(inputData, vl);
      inputData = __riscv_vxor_vv_u64m1_tu(inputData, inputData, crcAcc, 2); // vl=2 since we only want to XOR the 128-bit crcAcc 

      // Actual multiplication
      // Note: since the constant multiplicands are 32-bit wide, the upper 32-bit of vclmulh results are always 0
      vuint64m1_t multResLo = __riscv_vclmul_vv_u64m1(inputData, extRedCstVector, vl);
      vuint64m1_t multResHi = __riscv_vclmulh_vv_u64m1(inputData, extRedCstVector, vl);
      crcAccLo = __riscv_vredxor_vs_u64m1_u64m1(multResLo, zeroVecU64M1, vl);
      crcAccHi = __riscv_vredxor_vs_u64m1_u64m1(multResHi, zeroVecU64M1, vl);
      crcAcc = __riscv_vslideup_vx_u64m1(crcAccHi, crcAccLo, 1, 2);
  }
  // to finalize the CRC, we must:
  // 1. store the 128-bit of the folded value in memory and call the baseline scalar CRC on it
  // 2. finalize the CRC on the unprocessed bytes from p (should be 128-bit left)
  uint64_t crcAccBuffer[2] = {0, 0};
  // byte order reversing to ensure MS-byte is stored first (lowest address)
  crcAcc = __riscv_vrev8_v_u64m1(crcAcc, 2); 
  __riscv_vse64_v_u64m1(crcAccBuffer, crcAcc, 2);
  uint8_t* crcAccBufferU8 = (uint8_t*) crcAccBuffer;

  const uint32_t ethCRC32Poly = 0x04C11DB7;
  return crc32_be_generic(0, crcAccBufferU8, 16, ethCRC32Poly) ^ crc32_be_generic(0, p, len, ethCRC32Poly);
}

#if 0
uint32_t crcEth32_be_vector_opt_m2(uint32_t crc, unsigned char const *p, size_t len)
{
  uint32_t polynomial = 0x04C11DB7;
	int i;
  size_t avl = len / 8; // 8-byte per element

  // we take a pair of 64-bit elements (E, F) from the message and reduce
  // E.X^96 by multiplying  E.R where R=X^96[CRCPoly]
  // F.X^160 by multiplying F.S where S=X^160[CRCPoly]

   // {X^288, X^224, (X^160 mod P), (X^96 P)}
   const uint32_t redConstants[] = {0x569700e5, 0x75be46b7, 0xc5b9cd4c, 0xe8a45605};
   vuint32m1_t redConstantVector = __riscv_vle32_v_u32m1(redConstants, 2);
   vuint64m2_t extRedCstVector = __riscv_vzext_vf2_u64m2(redConstantVector, 2);

   vuint64m2_t crcAcc = __riscv_vmv_v_x_u64m2(0, 2);
   const vuint64m1_t zeroVecU64M1 = __riscv_vmv_v_x_u64m1(0, 2);
   const vuint64m2_t zeroVecU64M2 = __riscv_vmv_v_x_u64m2(0, 2);

   // hoisting vl
   size_t vl = __riscv_vsetvl_e64m2(-1);

  for (; avl >= 2*vl; avl -= vl, p += 8 * vl, len -= 8*vl) {
      // compute loop body vector length from application vector length avl
      vuint64m2_t inputData = __riscv_vle64_v_u64m2((uint64_t*) p, vl);
      // byte swapping the data to align their endianess with the CRC accumulator
      inputData = __riscv_vrev8_v_u64m2(inputData, vl);
      inputData = __riscv_vxor_vv_u64m2_tu(inputData, inputData, crcAcc, 2); // vl=2 since we only want to XOR the 128-bit crcAcc 

      // Actual multiplication
      // Note: since the constant multiplicands are 32-bit wide, the upper 32-bit of vclmulh results are always 0
      vuint64m2_t multResLo = __riscv_vclmul_vv_u64m2(inputData, extRedCstVector, vl);
      vuint64m2_t multResHi = __riscv_vclmulh_vv_u64m2(inputData, extRedCstVector, vl);
      vuint64m1_t crcAccLo = __riscv_vredxor_vs_u64m2_u64m1(multResLo, zeroVecU64M1, vl);
      vuint64m1_t crcAccHi = __riscv_vredxor_vs_u64m2_u64m1(multResHi, zeroVecU64M1, vl);
      crcAcc = __riscv_vslideup_vx_u64m2(crcAccHi, crcAccLo, 1, 2);
  }
  // to finalize the CRC, we must:
  // 1. store the 128-bit of the folded value in memory and call the baseline scalar CRC on it
  // 2. finalize the CRC on the unprocessed bytes from p (should be 128-bit left)
  uint64_t crcAccBuffer[2] = {0, 0};
  // byte order reversing to ensure MS-byte is stored first (lowest address)
  crcAcc = __riscv_vrev8_v_u64m1(crcAcc, 2); 
  __riscv_vse64_v_u64m1(crcAccBuffer, crcAcc, 2);
  uint8_t* crcAccBufferU8 = (uint8_t*) crcAccBuffer;

  const uint32_t ethCRC32Poly = 0x04C11DB7;
  return crc32_be_generic(0, crcAccBufferU8, 16, ethCRC32Poly) ^ crc32_be_generic(0, p, len, ethCRC32Poly);
}
#endif