#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <riscv_vector.h>

/** declaration of generic crc32_le_generic function, used to complete
 *  the vector implementations (epilog).
*/
uint32_t crc32_le_generic(uint32_t crc, unsigned char const *p,
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
uint32_t crcEth32_le_vector(uint32_t crc, unsigned char const *p, size_t len)
{
  uint32_t polynomial = 0x04C11DB7;
	int i;
  size_t avl = len / 4; // 4-byte per element 

  // we take a pair of 32-bit elements (E, F) from the message and reduce
  // E.X^32 by multiplying E.R where R=X^32[CRCPoly]
  // F.X^64 by multiplying F.S where S=X^64[CRCPoly]

  // {(X^96 mod P), (X^64 mod P)} = {0xe8a45605, 0xf200aa66};
  // {(X^64 mod P), (X^32 mod P)}
  // const uint32_t redConstants[] = {0x177b1443, 0x3d6029b0};
  const uint32_t redConstants[] = {
    /* CRC32(X^64) */ 0xf200aa66,
    /* CRC32(X^32) */ 0x490d678d,
  };
  vuint32mf2_t redConstantVector = __riscv_vle32_v_u32mf2(redConstants, 2);
  vuint64m1_t extRedCstVector = __riscv_vzext_vf2_u64m1(redConstantVector, 2);

  vuint64m1_t crcAcc = __riscv_vmv_v_x_u64m1(0, 2);
  const vuint64m1_t zeroVecU64M1 = __riscv_vmv_v_x_u64m1(0, 2);

  // hoisting vl value setting outside of loop (require test to check
  // that avl is always larger than 2*vl: one time for loop body and one
  // time for the epilog which is not vectorized).
  // size_t vl = __riscv_vsetvl_e32mf2(-1);

  for (size_t vl = __riscv_vsetvl_e32mf2(-1); avl >= 2 * vl; avl -= vl, p += 4 * vl, len -= 4*vl) {
      // compute loop body vector length from application vector length avl
      // printf("vector loop body, vl=%lu, avl=%lu, len=%lu\n", vl, avl, len);
      vuint32mf2_t inputData = __riscv_vle32_v_u32mf2((uint32_t*) p, vl);
      vuint32m1_t inputDataM1 = __riscv_vlmul_ext_v_u32mf2_u32m1(inputData);
      // printf("pre brev inputData=%"PRIx64"\n", __riscv_vmv_x_s_u64m1_u64(__riscv_vreinterpret_v_u32m1_u64m1(inputDataM1)));
      // bit swapping each 32-bit element
#if HAS_ZVBB_SUPPORT
      crcAcc = __riscv_vbrev_v_u64m1(crcAcc, 1); 
#else // no HAS_ZVBB_SUPPORT
      crcAcc = __riscv_vbrev8_v_u64m1(crcAcc, 1); 
      crcAcc = __riscv_vrev8_v_u64m1(crcAcc, 1); 
#endif
      //printf("post brev crcAccU64=%"PRIx64"\n", __riscv_vmv_x_s_u64m1_u64(crcAcc));
      vuint32m1_t crcAccU32 = __riscv_vreinterpret_v_u64m1_u32m1 (crcAcc);
      //printf("post brev crcAccU32=%"PRIx32"\n", __riscv_vmv_x_s_u32m1_u32(crcAccU32));
      //printf("post brev crcAccU32=%"PRIx32"\n", __riscv_vmv_x_s_u32m1_u32(__riscv_vslidedown_vx_u32m1(crcAccU32, 1, 2)));
      vuint32mf2_t crcAccU32mf2 = __riscv_vlmul_trunc_v_u32m1_u32mf2(crcAccU32);
      inputData = __riscv_vxor_vv_u32mf2(inputData, crcAccU32mf2, vl);
#if HAS_ZVBB_SUPPORT
      inputData = __riscv_vbrev_v_u32mf2(inputData, vl);
#else // no HAS_ZVBB_SUPPORT
      inputData = __riscv_vbrev8_v_u32mf2(inputData, vl);
      inputData = __riscv_vrev8_v_u32mf2(inputData, vl);
#endif
      // inputDataM1 = __riscv_vlmul_ext_v_u32mf2_u32m1(inputData);
      // printf("post xor+brev inputData=%"PRIx64"\n", __riscv_vmv_x_s_u64m1_u64(__riscv_vreinterpret_v_u32m1_u64m1(inputDataM1)));

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
  // bit order reversing to ensure MS-byte is stored first (lowest address)
#if HAS_ZVBB_SUPPORT
      crcAcc = __riscv_vbrev_v_u64m1(crcAcc, 1); 
#else // no HAS_ZVBB_SUPPORT
      crcAcc = __riscv_vbrev8_v_u64m1(crcAcc, 1); 
      crcAcc = __riscv_vrev8_v_u64m1(crcAcc, 1); 
#endif
  __riscv_vse64_v_u64m1(crcAccBuffer, crcAcc, 1);
  uint8_t* crcAccBufferU8 = (uint8_t*) crcAccBuffer;

  // const uint32_t ethCRC32Poly = 0x04C11DB7;
  const uint32_t ethCRC32PolyInv = 0xedb88320;
  return crc32_le_generic(0, crcAccBufferU8, 8, ethCRC32PolyInv) ^ crc32_le_generic(0, p, len, ethCRC32PolyInv);
}

