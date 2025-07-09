import math
import random
import argparse


def testCaseDescriptor(lhsSign, rhsSign, dividend, divisor):
     desc = f"({'neg' if lhsSign else 'pos'}-{'neg' if rhsSign else 'pos'}) divisor={hex(divisor)}, delta_ldc={leadingDigitPos - divisorLDC}"
     case = f"{{.v={{{hex(dividend)}ull, {hex(divisor)}ull}}, .label=\"{desc}\"}},"
     return case

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate test cases for 2-operand integer operations.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()
    
    random.seed(args.seed)
    print("generated_data_t data_2op_int[] = {")
    for lhsSign in [False, True]:
        for rhsSign in [False, True]:
            for divisor in [3, 0xcafe, 0xdeadbeef]:
                divisorLDC = int(math.ceil(math.log2(divisor)))
                divisor = -divisor if rhsSign else divisor
                for leadingDigitPos in range(divisorLDC, 64 if lhsSign else 65):
                    dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
                    dividend = -dividend if lhsSign else dividend
                    case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor)
                    print(f"  {case}")
    lhsSign, rhsSign = False, False
    divisor = 0x200 # 2^9
    divisorLDC = int(math.ceil(math.log2(divisor)))
    divisor = -divisor if rhsSign else divisor
    for leadingDigitPos in range(divisorLDC, 64 if lhsSign else 65):
        dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
        dividend = -dividend if lhsSign else dividend
        case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor)
        print(f"  {case}")
    divisor = 0x40 # 2^6
    divisorLDC = int(math.ceil(math.log2(divisor)))
    divisor = -divisor if rhsSign else divisor
    for leadingDigitPos in range(divisorLDC, 64 if lhsSign else 65):
        dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
        dividend = -dividend if lhsSign else dividend
        case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor)
        print(f"  {case}")
    print("};")