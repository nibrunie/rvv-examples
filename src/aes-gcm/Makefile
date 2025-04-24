
CRC_BE_BITS ?= 1

CLANG ?= clang
PK ?= /home/riscv_srcs/riscv-gnu-toolchain/pk/pk

INCLUDE ?= "-I/opt/riscv/riscv64-unknown-elf/include/"

SPIKE ?= spike

build/aes-gcm-test.o:
	mkdir -p build
	$(CLANG) -c --target=riscv64  -menable-experimental-extensions $(INCLUDE) -march=rv64gcv_zvbc1p0_zvbb1p0_zvbc32e1p0_zvkg_zvkgs1p0_zvkned aes-gcm-test.c $(CFLAGS) -o $@

build/zvkgs.o: zvkgs.c
	mkdir -p build
	$(CLANG) -c --target=riscv64  -menable-experimental-extensions $(INCLUDE) -march=rv64gcv_zvbc1p0_zvbb1p0_zvbc32e1p0_zvkg_zvkgs1p0_zvkned aes-gcm-test.c $(CFLAGS) -o $@
	
aes_gcm_test: build/aes-gcm-test.o zvkg.s zvb-ghash.s build/zvkgs.o log.c vlen-bits.s zvkned.s
	mkdir -p build
	riscv64-unknown-elf-gcc -DCRC_BE_BITS=$(CRC_BE_BITS) -I../utils -march=rv64gcv_zvkg_zvbc_zvbb_zvkned $^ -o $@ $(CFLAGS)

clean:
	rm *.o build/*.o aes_gcm_test

sim: aes_gcm_test
	$(SPIKE) --isa=rv64gcv_zvbc_zicntr_zihpm_zvbb_zvkned_zvkg_zvbc_zvbc32e_zvkgs $(PK) aes_gcm_test

.PHONY: sim