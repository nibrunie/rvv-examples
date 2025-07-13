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
            # one divisor per possible width
            for divisorBits in range(2, 65):
                for leadingDigitPos in range(2, 64 if lhsSign else 65):
                    for numCases in range(1, 2):
                        divisor = random.randrange(2**(divisorBits - 1), 2**divisorBits)
                        divisorLDC = int(math.ceil(math.log2(divisor)))
                        divisor = (2**64 - divisor) if rhsSign else divisor
                        # one dividend per possible width
                        dividend = (1 << (leadingDigitPos - 1)) | random.randrange(2** (leadingDigitPos - 1))
                        dividend = (2**64 - dividend) if lhsSign else dividend
                        case = testCaseDescriptor(lhsSign, rhsSign, dividend, divisor, leadingDigitPos, divisorLDC)
                        print(f"  {case}")
    print("};")
