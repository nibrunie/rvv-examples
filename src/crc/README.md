# Implementing CRC with RISC-V Vector

https://github.com/nibrunie/rvv-examples/tree/master/src/crc

This repository contains a few toy examples of CRC (Cyclic Redundancy Check) implementations, including implementations relying on RISC-V Vector official Zvbc extension and an implementation relying on an experimental new extension: Zvbc32e.

Those examples are described in a blog post on https://fprox.substack.com.

# Using Docker as development environment

## Building the image

The examples in this directory require an experimental RISC-V.
This toolchain can be build inside a Docker container using one of the docker file
present in this repository: https://github.com/nibrunie/rvv-examples/tree/master/riscv-toolchain-hub.Dockerfile

```
# from the top directory of rvv-examples repository
docker build -t riscv:rvv-examples-crc -f riscv-toolchain-hub.Dockerfile . 
```


## Running the image in a container

```
# from the top directory of rvv-examples repository
docker run  -ti --mount type=bind,source="$(pwd)"/,target=/home/app riscv:rvv-examples-crc
```

# Building and executing CRC benchmarks

From within the docker container previous built, you can run the CRC benchmarks easily:

```
cd src/crc
make clean; make SPIKE=/home/riscv_srcs/spike-vka/build/spike PK=/home/riscv_srcs/riscv-gnu-toolchain/build-pk64/pk sim CFLAGS=" -O3 -DCOUNT_INSTRET -DVERBOSE"
```

# Author(s)

Nicolas Brunie
