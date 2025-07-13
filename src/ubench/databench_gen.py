import math
import random
import argparse


def testCaseDescriptor(lhsSign, rhsSign, dividend, divisor, dividendLDC, divisorLDC):
     desc = f"({'neg' if lhsSign else 'pos'}-{'neg' if rhsSign else 'pos'}), {hex(dividend)}, {hex(divisor)}, {dividendLDC - divisorLDC}, {dividendLDC}, {divisorLDC}, {hex(dividend // divisor)}"
     case = f"{{.v={{{hex(dividend)}ull, {hex(divisor)}ull}}, .label=\"{desc}\"}},"
     return case

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate test cases for 2-operand integer operations.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()
    
    random.seed(args.seed)
    print("generated_data_t data_2op_int[] = {")
    for lhsSign in [False]: # , True]:
        for rhsSign in [False]: # , True]:
            for divisor in [3, 0xcafe, 0x3cafe, 0x7eadbeef, 0x7eadcafebeef] + [random.randrange(3, 2**random.randrange(2, 65)) for i in range(1000)]:
                divisorLDC = int(math.ceil(math.log2(divisor)))
                divisor = (2**64 - divisor) if rhsSign else divisor
                leadingDigitPos = random.randrange(divisorLDC, 64 if lhsSign else 65)
                # for leadingDigitPos in range(divisorLDC, 64 if lhsSign else 65):
                if True:
                    dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
                    dividend = (2**64 - dividend) if lhsSign else dividend
                    case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor, leadingDigitPos, divisorLDC)
                    print(f"  {case}")
    lhsSign, rhsSign = False, False
    divisor = 0x200 # 2^9
    divisorLDC = int(math.ceil(math.log2(divisor)))
    divisor = -divisor if rhsSign else divisor
    for leadingDigitPos in range(divisorLDC, 64 if lhsSign else 65):
        dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
        dividend = -dividend if lhsSign else dividend
        case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor, leadingDigitPos, divisorLDC)
        print(f"  {case}")
    divisor = 0x40 # 2^6
    divisorLDC = int(math.ceil(math.log2(divisor)))
    divisor = -divisor if rhsSign else divisor
    for leadingDigitPos in range(divisorLDC, 64 if lhsSign else 65):
        dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
        dividend = -dividend if lhsSign else dividend
        case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor, leadingDigitPos, divisorLDC)
        print(f"  {case}")
    for divisor in [0x40, 3, 0xcafe, 0x17, 0x1eadbeef, 0x3eadbeef, 0x5eadbeeef, 0xdeadbeef, 0x1deadbeef]:
        divisorLDC = int(math.ceil(math.log2(divisor)))
        divisor = -divisor if rhsSign else divisor
        for leadingDigitPos in range(divisorLDC, 63 if lhsSign else 64):
            # special case: dividend is a power of two times the divisor                                                                              
            scale = 2**(leadingDigitPos  - divisorLDC)
            dividend = divisor * scale
            dividend = -dividend if lhsSign else dividend
            case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor, leadingDigitPos, divisorLDC)
            print(f"  {case}")
            # random case
            #dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
            #dividend = (2**64 - dividend) if lhsSign else dividend
            #case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor, leadingDigitPos, divisorLDC)
            #print(f"  {case}")
    print("};")
