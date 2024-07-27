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

/** crc32_be_vector_zvbc32e() - Calculate bitwise little-endian Ethernet AUTODIN II
 *			CRC32/CRC32C using Zvbc32e
 * @crc: seed value for computation.  ~0 for Ethernet, sometimes 0 for other
 *	 uses, or the previous crc32/crc32c value if computing incrementally.
 * @p: pointer to buffer over which CRC32/CRC32C is run
 * @len: length of buffer @p
 * @redCsts: reduction consts {X^160, X^96, ...}
 *
 */
uint32_t crc_be_vector_zvbc32e(uint32_t crc, unsigned char const *p, size_t len, const uint32_t redCsts[])
{
	int i;
  size_t avl = len / 4; // 4-byte per element


  // we take 4 32-bit elements (E, F, G, H) from the message and reduce
  // them 
  // E.X^64  by multiplying E.R where R=X^64 [CRCPoly]
  // F.X^96  by multiplying F.S where S=X^96 [CRCPoly]
  // G.X^128 by multiplying F.T where S=X^128[CRCPoly]
  // H.X^160 by multiplying F.U where S=X^160[CRCPoly]
  // we expect {(X^160 mod P), (X^128 mod P), (X^96 mod P), (X^64 mod P)}
  // to be stored in redCsts table

   vuint32m1_t redConstantVector = __riscv_vle32_v_u32m1(redCsts, 4);
   for (int i = 0; i < 4; ++i) {
    printf("redCsts[%d] = %"PRIx32"\n", i, redCsts[i]);
   }

   vuint32m1_t crcAccLo = __riscv_vmv_v_x_u32m1(0, 4);
   vuint32m1_t crcAccHi = __riscv_vmv_v_x_u32m1(0, 4);
   vuint32m1_t crcAcc   = __riscv_vmv_v_x_u32m1(0, 4);
   const vuint32m1_t zeroVecU32M1 = __riscv_vmv_v_x_u32m1(0, 4);

   // hoisting setting vl (assuming it won't change in loop)
   size_t vl = __riscv_vsetvl_e32m1(-1);

  for (; avl >= 2*vl; avl -= vl, p += 4 * vl, len -= 4*vl) {
      printf("loop avl=%lu.\n", avl);
      // compute loop body vector length from application vector length avl
      assert(vl == 4);
      vuint32m1_t inputData = __riscv_vle32_v_u32m1((uint32_t*) p, vl);
      // byte swapping the data to align their endianess with the CRC accumulator
      inputData = __riscv_vrev8_v_u32m1(inputData, vl);
      inputData = __riscv_vxor_vv_u32m1_tu(inputData, inputData, crcAcc, 2); // vl=2 since we only want to XOR the 64-bit crcAcc 

      // Actual multiplication
      vuint32m1_t multResLo; // = __riscv_vclmul_vv_u32m1(inputData, extRedCstVector, vl);
      vuint32m1_t multResHi; // = __riscv_vclmulh_vv_u32m1(inputData, extRedCstVector, vl);
      asm volatile ("vclmul.vv %0, %1, %2" : "=vr"(multResLo) : "vr"(inputData), "vr"(redConstantVector));
      asm volatile ("vclmulh.vv %0, %1, %2" : "=vr"(multResHi) : "vr"(inputData), "vr"(redConstantVector));

      crcAccLo = __riscv_vredxor_vs_u32m1_u32m1(multResLo, zeroVecU32M1, vl);
      crcAccHi = __riscv_vredxor_vs_u32m1_u32m1(multResHi, zeroVecU32M1, vl);
      crcAcc = __riscv_vslideup_vx_u32m1(crcAccHi, crcAccLo, 1, 2);
  }
  printf("final avl=%lu.\n", avl);
  // to finalize the CRC, we must:
  // 1. store the 64-bit of the folded value in memory and call the baseline scalar CRC on it
  // 2. finalize the CRC on the unprocessed bytes from p (should be 128-bit left)
  uint32_t crcAccBuffer[4] = {0, 0, 0, 0};
  // byte order reversing to ensure MS-byte is stored first (lowest address)
  crcAcc = __riscv_vrev8_v_u32m1(crcAcc, 2); 
  __riscv_vse32_v_u32m1(crcAccBuffer, crcAcc, 2);
  uint8_t* crcAccBufferU8 = (uint8_t*) crcAccBuffer;

  const uint32_t ethCRC32Poly = 0x04C11DB7;
  const uint32_t finalSeed = crc32_be_generic(0, crcAccBufferU8, 16, ethCRC32Poly);
  printf("finalSeed=%"PRIx32", len=%lu\n", finalSeed, len);
  return  finalSeed ^ crc32_be_generic(0, p, len, ethCRC32Poly);
}