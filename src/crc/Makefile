
CRC_BE_BITS ?= 1

CLANG = clang
PK ?= /home/riscv_srcs/riscv-gnu-toolchain/pk/pk

INCLUDE ?= "-I/opt/riscv/riscv64-unknown-elf/include/"

build/vector_crc_be.o:
	mkdir -p build
	$(CLANG) -c --target=riscv64  -menable-experimental-extensions $(INCLUDE) -march=rv64gcv_zvbc1p0_zvbb1p0 vector_crc_be.c $(CFLAGS) -o $@
build/vector_crc_le.o:
	mkdir -p build
	$(CLANG) -c --target=riscv64  -menable-experimental-extensions $(INCLUDE) -march=rv64gcv_zvbc1p0_zvbb1p0 vector_crc_le.c $(CFLAGS) -o $@

bench_crc: build/vector_crc_be.o build/vector_crc_le.o
	riscv64-unknown-elf-gcc -DCRC_BE_BITS=$(CRC_BE_BITS) -march=rv64gcv crc32.c crc_bench.c $^ -o $@ $(CFLAGS)

clean:
	rm *.o build/*.o bench_crc

sim: bench_crc
	spike --isa=rv64gcv_zvbc_zicntr_zihpm_zvbb $(PK) bench_crc

.PHONY: sim_crc