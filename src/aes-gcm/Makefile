
CRC_BE_BITS ?= 1

CLANG ?= clang
PK ?= /home/riscv_srcs/riscv-gnu-toolchain/pk/pk

INCLUDE ?= "-I/opt/riscv/riscv64-unknown-elf/include/"

SPIKE ?= spike

build/aes-gcm-test.o: aes-gcm-test.c
	mkdir -p build
	$(CLANG) -c --target=riscv64  -menable-experimental-extensions $(INCLUDE) -march=rv64gcv_zvbc1p0_zvbb1p0_zvbc32e1p0_zvkg_zvkgs1p0_zvkned $^ $(CFLAGS) -o $@

build/zvkgs.o: zvkgs.c
	mkdir -p build
	$(CLANG) -c --target=riscv64  -menable-experimental-extensions $(INCLUDE) -march=rv64gcv_zvbc1p0_zvbb1p0_zvbc32e1p0_zvkg_zvkgs1p0_zvkned $^ $(CFLAGS) -o $@
	
aes_gcm_test: build/aes-gcm-test.o zvkg.s zvb-ghash.s build/zvkgs.o log.c vlen-bits.s zvkned.s
	mkdir -p build
	riscv64-unknown-elf-gcc -I../utils -march=rv64gcv_zvkg_zvbc_zvbb_zvkned $^ -o $@ $(CFLAGS)

build/zvbc32e.o: zvbc32e.s
	mkdir -p build
	$(CLANG) -c --target=riscv64  -menable-experimental-extensions $(INCLUDE) -march=rv64gcv_zvbc1p0_zvbb1p0_zvbc32e1p0_zvkg_zvkgs1p0_zvkned $^ $(CFLAGS) -o $@

zvbc32e_test: build/zvbc32e.o zvb-ghash.s zvbc32e-test.c log.c vlen-bits.s zvkned.s
	mkdir -p build
	riscv64-unknown-elf-gcc -I../utils -march=rv64gcv_zvkg_zvbc_zvbb_zvkned $^ -o $@ $(CFLAGS)

clean:
	rm *.o build/*.o aes_gcm_test

sim_aes_gcm_zvkgs: aes_gcm_test
	$(SPIKE) --isa=rv64gcv_zvbc_zicntr_zihpm_zvbb_zvkned_zvkg_zvbc_zvbc32e_zvkgs $(PK) aes_gcm_test

sim_zvbc32e_test: zvbc32e_test
	$(SPIKE) --isa=rv64gcv_zvbc_zicntr_zihpm_zvbb_zvkned_zvkg_zvbc_zvbc32e_zvkgs $(PK) zvbc32e_test

.PHONY: sim