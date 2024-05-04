# Makefile for RISC-V In a Nutshell course (RVVIAN) examples

# RISC-V C Compiler
# available options (the version used must support RVV intrinsics)
# clang/llvm
RISCVCLANG=clang  --target=riscv64
# RISCVCC=clang  --target=riscv64
# GNU Compiler Collection (GCC)
RISCVCC=riscv64-unknown-elf-gcc 

EXTRA_CFLAGS?=-DCOUNT_INSTRET

# architectural parameters for the simulation
# width of vector registers (VLEN)
VLEN?=128

# path to proxy-kernel (pk)
PK_PATH=/opt/riscv/riscv64-unknown-elf/bin/pk64 

# SIMULATOR
# Available options in the Docker (uncomment one)
SIMULATOR=spike --isa=rv64gcv_zicntr_zihpm --varch=vlen:$(VLEN),elen:64 $(PK_PATH)
# SIMULATOR=qemu-riscv64 -cpu rv64,v=on,vext_spec=v1.0,vlen=128,rvv_ta_all_1s=on

INCLUDE_DIR ?= /opt/riscv/riscv64-unknown-elf/include/

poly_mult_baseline.o: poly_mult_baseline.c
	$(RISCVCLANG) $(EXTRA_CFLAGS) -I./ -I$(INCLUDE_DIR)  -O2 -march=rv64gcv -c -o $@ $^

poly_mult_rvv.o: poly_mult_rvv.c
	$(RISCVCLANG) $(EXTRA_CFLAGS) -I./ -I$(INCLUDE_DIR)  -O2 -march=rv64gcv -c -o $@ $^

bench_poly_mult: bench_poly_mult.c poly_mult_baseline.o poly_mult_rvv.o
	 $(RISCVCC) $(EXTRA_CFLAGS) -I./ -O2 -march=rv64gcv $^ -lm -o $@

basic_poly_test: basic_poly_test.c poly_mult_baseline.o
	 $(RISCVCC) $(EXTRA_CFLAGS) -I./ -O2 -march=rv64gcv $^ -lm -o $@

sim_bench_poly_mult: bench_poly_mult
	$(SIMULATOR) $^

sim_basic_poly_test: basic_poly_test
	$(SIMULATOR) $^

clean:
	rm -f *.o bench_poly_mult  basic_poly_test

.PHONY: sim_bench_poly_mult sim_basic_poly_test clean
