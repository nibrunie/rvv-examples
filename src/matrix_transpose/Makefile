# Makefile for RISC-V In a Nutshell course (RVVIAN) examples

# RISC-V C Compiler
# available options (the version used must support RVV intrinsics)
# clang/llvm
RISCVCLANG=clang  --target=riscv64
# RISCVCC=clang  --target=riscv64
# GNU Compiler Collection (GCC)
RISCVCC=riscv64-unknown-elf-gcc 


# architectural parameters for the simulation
# width of vector registers (VLEN)
VLEN?=128

# path to proxy-kernel (pk)
PK_PATH=/opt/riscv/riscv64-unknown-elf/bin/pk64 

# SIMULATOR
# Available options in the Docker (uncomment one)
SIMULATOR=spike --isa=rv64gcv_zicntr_zihpm --varch=vlen:$(VLEN),elen:64 $(PK_PATH)
# SIMULATOR=qemu-riscv64 -cpu rv64,v=on,vext_spec=v1.0,vlen=128,rvv_ta_all_1s=on

matrix_transpose_intrinsics.o: matrix_transpose_intrinsics.c
	$(RISCVCLANG) -I./ -O2 -march=rv64gcv -c -o $@ $^

bench_matrix_transpose: bench_matrix_transpose.c matrix_transpose_intrinsics.o
	 $(RISCVCC) -I./ -O2 -march=rv64gcv $^ -o $@

sim_bench_matrix_transpose: bench_matrix_transpose
	$(SIMULATOR) $^

clean:
	rm -f *.o bench_matrix_transpose 

.PHONY: sim_bench_matrix_transpose clean
