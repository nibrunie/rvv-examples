# Makefile for RISC-V In a Nutshell course (RVVIAN) examples

# RISC-V C Compiler
# available options (the version used must support RVV intrinsics)
# clang/llvm
RISCVCC=clang  --target=riscv64
# GNU Compiler Collection (GCC)
# RISCVCC=riscv64-unknown-elf-gcc 

# architectural parameters for the simulation
# width of vector registers (VLEN)
VLEN?=128

# path to proxy-kernel (pk)
PK_PATH=/opt/riscv/riscv64-unknown-elf/bin/pk 

# SIMULATOR
# Available options in the Docker (uncomment one)
SIMULATOR=spike --isa=rv64gcv_zicntr_zihpm --varch=vlen:$(VLEN),elen:64 $(PK_PATH)
# SIMULATOR=qemu-riscv64 -cpu rv64,v=on,vext_spec=v1.0,vlen=128,rvv_ta_all_1s=on

vector_add_intrinsics.o: vector_add_intrinsics.c
	$(RISCVCC) -O2 -march=rv64gcv -c -o $@ $^

bench_vector_add: bench_vector_add.c vector_add_intrinsics.o
	 $(RISCVCC) -O2 -march=rv64gcv $^ -o $@

sim_bench_vector_add: bench-0_vector_add
	$(SIMULATOR) $^

clean:
	rm -f *.o bench_vector_add 

.PHONY: sim_bench_vector_add clean
