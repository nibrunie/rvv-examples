# rvv-examples
Example of RISC-V Vector programming

## Organization

# Building the docker image

This project relies on Docker to reliably build a working RISC-V build and run environment.
You can build the image with the following command:

```
docker build  -t riscv:riscv-toolchain -f riscv-toolchain.Dockerfile . 
```

The built image will be stored with the tag `riscv:riscv-toolchain`

# Running example in a docker container

The easiest way to build and run the example is to run them in a docker container started on the image previously built.
By running the following command while at the top of the `./rvv-examples` directory you can bind the full directory content to the directory `/home/app` within the container. Any modification made within the container will be visible outside of it (and the opposite is also true).

```
docker run  -ti --mount type=bind,source="$(pwd)"/,target=/home/app/ riscv:riscv-toolchain
```
