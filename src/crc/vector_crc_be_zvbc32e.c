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
uint32_t crcEth32_be_vector_zvbc32e(uint32_t crc, unsigned char const *p, size_t len)
{
	int i;
  size_t avl = len / 4; // 4-byte per element
  const uint32_t ethCRC32Poly = 0x04C11DB7;

  // fast exit if length is not large enough to justify entering the vectorized loop
  // FIXME: The threshold could be tuned
  if (len <= 16) return crc32_be_generic(0, p, len, ethCRC32Poly);

  // we take 4 32-bit elements (E, F, G, H) from the message and reduce
  // them 
  // E.X^64  by multiplying E.R where R=X^64 [CRCPoly]
  // F.X^96  by multiplying F.S where S=X^96 [CRCPoly]
  // G.X^128 by multiplying F.T where S=X^128[CRCPoly]
  // H.X^160 by multiplying F.U where S=X^160[CRCPoly]
  // {(X^160 mod P), (X^128 mod P), (X^96 mod P), (X^64 mod P)}
  // are stored in redCsts table
  // Those constants were built as follow
  // {
  //   uint32_t ethCRC32Poly = 0x04C11DB7;
  //   uint8_t dbgMsg[17] = {0};
  //   dbgMsg[0] = 0x1;
  //   for (int i = 0; i < 4; ++i) redCsts[i] = crc32_be_generic(0, dbgMsg, 1 + 16 - 4 * i, ethCRC32Poly);
  // }
  const uint32_t redCsts[4] = {
    0x17d3315d,
    0xe8a45605,
    0xf200aa66,
    0x490d678d,
  };

   vuint32m1_t redConstantVector = __riscv_vle32_v_u32m1(redCsts, 4);

   vuint32m1_t crcAccLo = __riscv_vmv_v_x_u32m1(0, 4);
   vuint32m1_t crcAccHi = __riscv_vmv_v_x_u32m1(0, 4);
   vuint32m1_t crcAcc   = __riscv_vmv_v_x_u32m1(0, 4);
   const vuint32m1_t zeroVecU32M1 = __riscv_vmv_v_x_u32m1(0, 4);

   // hoisting setting vl (assuming it won't change in loop)
   size_t vl = __riscv_vsetvl_e32m1(-1);

#if 1
  asm volatile (
    "mv a2, %[bound]\n"
    "bgeu	%[p], a2, 2f\n" // skipping inner loop if length is not large enough
    "vmv.v.i v11, 0\n"
    "vsetivli	zero,4,e32,m1,tu,ma\n"
    "1:\n"
    "vle32.v	v13, (%[p])\n" // s1 -> p
    "vrev8.v	v13, v13\n"
    "vxor.vv	v13, v13, %[crcAcc]\n" // v10 -> crcAcc
    "vclmul.vv	v20, v13, %[redConstantVector]\n"
    "vclmulh.vv	v13, v13, %[redConstantVector]\n"
    "vredxor.vs	v11, v20, v8\n"
    "vredxor.vs	%[crcAcc], v13, v8\n"
    "vslideup.vi	%[crcAcc], v11, 1\n"
    "add	%[p], %[p], 16\n"
    "bltu	%[p], a2, 1b\n"
    "2:\n"
  : [p]"+r"(p), [len]"+r"(len), [avl]"+r"(avl), [crcAcc]"+vr"(crcAcc), [redConstantVector]"+vr"(redConstantVector)
  // : [bound]"r"(2*vl*4)
  : [bound]"r"(p + len - 16) // FIXME bound assume len is a multiple of 16
  : "v10", "v13", "v20", "v12", "v8", "v11", "a2"
  );
  len = 16; // remainder after end of loop (FIXME: assume original len was a multiple of 16)

#else
  for (; avl >= 2*vl; avl -= vl, p += 4 * vl, len -= 4*vl) {
      // printf("loop avl=%lu.\n", avl);
      // compute loop body vector length from application vector length avl
      assert(vl == 4);
      vuint32m1_t inputData = __riscv_vle32_v_u32m1((uint32_t*) p, vl);
      // byte swapping the data to align their endianess with the CRC accumulator
      inputData = __riscv_vrev8_v_u32m1(inputData, vl);
      inputData = __riscv_vxor_vv_u32m1_tu(inputData, inputData, crcAcc, 2); // vl=2 since we only want to XOR the 64-bit crcAcc 

      // Actual multiplication
      vuint32m1_t multResLo; // = __riscv_vclmul_vv_u32m1(inputData, extRedCstVector, vl);
      vuint32m1_t multResHi; // = __riscv_vclmulh_vv_u32m1(inputData, extRedCstVector, vl);
      asm volatile (
        "vsetivli x0, 4, e32, m1, tu, mu\n"
        "vclmul.vv v20, %[inputData], %[redConstantVector]\n"
        "vclmulh.vv %[multResultHi], %[inputData], %[redConstantVector]\n"
        "vmv.v.v %[multResultLo], v20\n"
         : [multResultLo]"=vr"(multResLo), [multResultHi]"=vr"(multResHi) 
         : [inputData]"vr"(inputData), [redConstantVector]"vr"(redConstantVector)
         : "v20");
      // asm volatile ("vclmulh.vv %0, %1, %2" : "=vr"(multResHi) : "vr"(inputData), "vr"(redConstantVector));

      crcAccLo = __riscv_vredxor_vs_u32m1_u32m1_tu(crcAccLo, multResLo, zeroVecU32M1, vl);
      crcAccHi = __riscv_vredxor_vs_u32m1_u32m1_tu(crcAccHi, multResHi, zeroVecU32M1, vl);
      crcAcc = __riscv_vslideup_vx_u32m1(crcAccHi, crcAccLo, 1, 2);
  }
#endif
  // to finalize the CRC, we must:
  // 1. store the 64-bit of the folded value in memory and call the baseline scalar CRC on it
  // 2. finalize the CRC on the unprocessed bytes from p (should be 128-bit left)
  uint32_t crcAccBuffer[4] = {0, 0, 0, 0};
  if (len >= 16)
  // byte order reversing to ensure MS-byte is stored first (lowest address)
  crcAcc = __riscv_vrev8_v_u32m1(crcAcc, 2); 
  __riscv_vse32_v_u32m1(crcAccBuffer, crcAcc, 2);
  uint8_t* crcAccBufferU8 = (uint8_t*) crcAccBuffer;

  const uint32_t finalSeed = crc32_be_generic(0, crcAccBufferU8, 16, ethCRC32Poly);
  return  finalSeed ^ crc32_be_generic(0, p, len, ethCRC32Poly);
}