# Makefile for RISC-V In a Nutshell course (RVVIAN) examples

# RISC-V C Compiler
# available options (the version used must support RVV intrinsics)
# clang/llvm
CLANG_EXTRA_EXTS=_zvbb1p0
RISCVCC=clang  --target=riscv64 -I/opt/riscv/riscv64-unknown-elf/include/
# GNU Compiler Collection (GCC)
# RISCVCC=riscv64-unknown-elf-gcc 

RISCVGNUCC=riscv64-unknown-elf-gcc 

# architectural parameters for the simulation
# width of vector registers (VLEN)
VLEN?=128

# path to proxy-kernel (pk)
PK_PATH=/opt/riscv/riscv64-unknown-elf/bin/pk64 

# SIMULATOR
# Available options in the Docker (uncomment one)
SIMULATOR=spike --isa=rv64gcv_zicntr_zihpm_zvbb --varch=vlen:$(VLEN),elen:64 $(PK_PATH)
# SIMULATOR=qemu-riscv64 -perfmap -cpu rv64,zba=true,v=on,vext_spec=v1.0,vlen=128,rvv_ta_all_1s=on

CFLAGS ?= -O2

keccak.o: keccak-vector.c # Keccak-more-compact.c
	 $(RISCVCC) $(CFLAGS) -c -menable-experimental-extensions -march=rv64gcv$(CLANG_EXTRA_EXTS) $^ -o $@

bench_keccak: Keccak-compact-test.c keccak.o genKAT.c
	 $(RISCVGNUCC) $(CFLAGS) -march=rv64gcv $^ -o $@

sim_bench_keccak: bench_keccak
	$(SIMULATOR) $^

clean:
	rm -f *.o bench_keccak 

.PHONY: sim_bench_keccak clean
