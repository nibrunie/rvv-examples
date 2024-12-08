# Makefile for RISC-V In a Nutshell course (RVVIAN) examples

# RISC-V C Compiler
# available options (the version used must support RVV intrinsics)
# clang/llvm
RISCVCLANG=clang  --target=riscv64
# RISCVCC=clang  --target=riscv64
# GNU Compiler Collection (GCC)
RISCVCC=riscv64-unknown-elf-gcc 
OBJDUMP=riscv64-unknown-elf-objdump

EXTRA_CFLAGS?=-DCOUNT_INSTRET -O2

# architectural parameters for the simulation
# width of vector registers (VLEN)
VLEN?=128

# path to proxy-kernel (pk)
PK_PATH=/opt/riscv/riscv64-unknown-elf/bin/pk64 

BUILD_ARCH=rv64gcv_zvbb1p0 -menable-experimental-extensions
SIM_ARCH=rv64gcv_zvbb

# SIMULATOR
# Available options in the Docker (uncomment one)
SIMULATOR=spike --isa=$(SIM_ARCH)_zicntr_zihpm --varch=vlen:$(VLEN),elen:64 $(PK_PATH)
# SIMULATOR=qemu-riscv64 -cpu rv64,v=on,vext_spec=v1.0,vlen=128,rvv_ta_all_1s=on

INCLUDE_DIR ?= /opt/riscv/riscv64-unknown-elf/include/
BUILD_DIR ?= $(PWD)/build/

$(BUILD_DIR)/poly_mult_baseline.o: poly_mult_baseline.c
	mkdir -p $(BUILD_DIR)
	$(RISCVCLANG) $(EXTRA_CFLAGS) -I./ -I$(INCLUDE_DIR) -march=$(BUILD_ARCH) -c -o $@ $^

$(BUILD_DIR)/poly_mult_rvv.o: poly_mult_rvv.c
	mkdir -p $(BUILD_DIR)
	$(RISCVCLANG) $(EXTRA_CFLAGS) -I./ -I$(INCLUDE_DIR) -march=$(BUILD_ARCH) -c -o $@ $^

$(BUILD_DIR)/ntt_scalar.o: ntt_scalar.c
	mkdir -p $(BUILD_DIR)
	$(RISCVCLANG) $(EXTRA_CFLAGS) -I./ -I$(INCLUDE_DIR) -march=$(BUILD_ARCH) -c -o $@ $^

$(BUILD_DIR)/poly_mult_ntt.o: poly_mult_ntt.c
	mkdir -p $(BUILD_DIR)
	$(RISCVCLANG) $(EXTRA_CFLAGS) -I./ -I$(INCLUDE_DIR) -march=$(BUILD_ARCH) -c -o $@ $^

$(BUILD_DIR)/bench_poly_mult: bench_poly_mult.c $(BUILD_DIR)/poly_mult_baseline.o $(BUILD_DIR)/poly_mult_rvv.o $(BUILD_DIR)/ntt_scalar.o $(BUILD_DIR)/poly_mult_ntt.o
	mkdir -p $(BUILD_DIR)
	$(RISCVCC) $(EXTRA_CFLAGS) -I./ -march=rv64gcv $^ -lm -o $@

$(BUILD_DIR)/basic_poly_test: basic_poly_test.c $(BUILD_DIR)/poly_mult_baseline.o $(BUILD_DIR)/ntt_scalar.o
	mkdir -p $(BUILD_DIR)
	$(RISCVCC) $(EXTRA_CFLAGS) -I./ -march=rv64gcv $^ -lm -o $@

sim_bench_poly_mult: $(BUILD_DIR)/bench_poly_mult
	$(SIMULATOR) $^

sim_basic_poly_test: $(BUILD_DIR)/basic_poly_test
	$(SIMULATOR) $^

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/bench_poly_mult $(BUILD_DIR)/basic_poly_test

.PHONY: sim_bench_poly_mult sim_basic_poly_test clean
