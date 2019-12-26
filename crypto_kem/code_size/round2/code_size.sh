#!/bin/bash

#kyber files
riscv64-unknown-elf-gcc -O3 -fno-common -o verify.o -c verify.c
riscv64-unknown-elf-gcc -O3 -fno-common -o indcpa.o -c indcpa.c
riscv64-unknown-elf-gcc -O3 -fno-common -o kem.o -c kem.c
riscv64-unknown-elf-gcc -O3 -fno-common -o nttc.o -c ntt.c
riscv64-unknown-elf-gcc -O3 -fno-common -o poly.o -c poly.c
riscv64-unknown-elf-gcc -O3 -fno-common -o polyvec.o -c polyvec.c
riscv64-unknown-elf-gcc -O3 -fno-common -o reduce.o -c reduce.c
riscv64-unknown-elf-gcc -O3 -fno-common -o symmetric-fips202.o -c symmetric-fips202.c
riscv64-unknown-elf-gcc -O3 -fno-common -o cbd.o -c cbd.c

#asm
riscv64-unknown-elf-gcc -O3 -fno-common -o ntt.o -c ntt_2.S
riscv64-unknown-elf-gcc -O3 -fno-common -o invntt.o -c invntt_2.S

#common files
riscv64-unknown-elf-gcc -O3 -fno-common -o fips202.o -c fips202.c
riscv64-unknown-elf-gcc -O3 -fno-common -o sha2.o -c sha2.c 
riscv64-unknown-elf-gcc -O3 -fno-common -o keccakf1600.o -c keccakf1600.c 
riscv64-unknown-elf-gcc -O3 -fno-common -o randombytes.o -c randombytes.c 


#size
riscv64-unknown-elf-ar -crs arch.a *.o
#riscv64-unknown-elf-nm arch.a 
echo "\n\n  KYBER SIZE ROUND 2:"
riscv64-unknown-elf-size -t arch.a
echo "\n\n"
#clean
rm *.o
rm *.a
