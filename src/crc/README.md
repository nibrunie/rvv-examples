# rvvian

RISC-V Vector In a Nutshell (examples for the course)

https://github.com/RVVIAN/rvvian


# Using Docker as development environment

## Requirements


## Building the image

```
cd DockerImage
docker build -t riscv:rvvian . 
```


## Running the image in a container

```
# from the top directory of rvvian repository
docker run  -ti --mount type=bind,source="$(pwd)"/Examples,target=/home/app riscv:rvvian
```

## Examples

### Example 0: vector addtion

#### Building

```
export PATH=/home/riscv_srcs/llvm-project/build/bin/:$PATH
clang -c --target=riscv64  -menable-experimental-extensions   -march=rv64gcv_zvbc1p0_zvbb1p0 src/vector_crc_intrinsics.c -O2
riscv64-unknown-elf-gcc -march=rv64gcv src/crc32.c src/crc_bench.c vector_crc_intrinsics.o -o bench_crc -O2
```

#### Executing with spike

Spike (https://github.com/riscv-software-src/riscv-isa-sim) is one of the main instruction set simulator for RISC-V.
It can be used with riscv proxy kernel, `pk`, (https://github.com/riscv-software-src/riscv-pk) to execute basic programs.

```
spike --isa=rv64gcv_zvbc_zicntr_zihpm /opt/riscv/riscv64-unknown-elf/bin/pk bench_crc
```

To generate assembly traces: add `-l` option.

#### Executing with QEMU

```
qemu-riscv64 -cpu rv64,v=on,vext_spec=v1.0,vlen=128,rvv_ta_all_1s=on ./bench_crc
```


# Author(s)
