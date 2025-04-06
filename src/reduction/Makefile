# Makefile for reduction benchmark

CLANG = clang
PK ?= /home/riscv_srcs/riscv-gnu-toolchain/pk/pk

INCLUDE ?= "-I/opt/riscv/riscv64-unknown-elf/include/"

SPIKE ?= spike


bench_reduction: reduction_bench.c
	riscv64-unknown-elf-gcc -O2 -I../utils -march=rv64gcv_zba $^ -o $@ $(CFLAGS)

clean:
	rm *.o build/*.o bench_reduction

sim: bench_reduction
	$(SPIKE) --isa=rv64gcv_zicntr_zihpm_zba $(PK) bench_reduction

.PHONY: sim