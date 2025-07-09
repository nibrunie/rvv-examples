import math
import random


def testCaseDescriptor(lhsSign, rhsSign, dividend, divisor):
     desc = f"({'neg' if lhsSign else 'pos'}-{'neg' if rhsSign else 'pos'}) divisor={hex(divisor)}, delta_ldc={leadingDigitPos - divisorLDC}"
     case = f".v={{{{{hex(dividend)}ull, {hex(divisor)}ull}}, .label=\"{desc}\"}},"
     return case

print("generated_data_t data_2op_int[] = {")
for lhsSign in [False, True]:
    for rhsSign in [False, True]:
        for divisor in [3, 0xcafe, 0xdeadbeef]:
            divisorLDC = int(math.ceil(math.log2(divisor)))
            for leadingDigitPos in range(divisorLDC, 65):
                dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
                case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor)
                print(f"  {case}")
    print("};")