# Kyber on RISC-V

**Author:** Denisa Greconici `<D.Greconici@cs.ru.nl>` \
This project is based on the [mupq/pqriscv](https://github.com/mupq/pqriscv) repository.\
The optimized Kyber implementation for RISC-V is in the folder:  `krypto_kem`


### Prerequisites
1. Make sure that you install the [RISC-V GNU toolchain](https://github.com/riscv/riscv-gnu-toolchain). 
2. Make sure your current jdk is 1.8 (java -version). It does not work with newer versions of jdk. If your jdk version is newer, run the following commands:
    ```sh
    $ sudo apt-get install openjdk-8-jre
    $ sudo update-alternatives –config java
    ```
    Choose java 8 from the list generated in the last command.
3. Install [SBT](https://www.scala-sbt.org/).
    For Ubuntu or Debian-based distribution, you can run the following commands:
    ```sh
    $ echo "deb https://dl.bintray.com/sbt/debian /" | sudo tee -a /etc/apt/sources.list.d/sbt.list
    $ curl -sL "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x2EE0EA64E40A89B84B2DF73499E82A75642AC823” | sudo apt-key add
    $ sudo apt-get update
    $ sudo apt-get install sbt
    ```
4. Install [Verilator](https://www.veripool.org/wiki/verilator).
5. Clone the thesis project:
    `$ git clone https://github.com/denigreco/Kyber_RISC_V_Thesis.git`
6. Clone the pqriscv-vexriscv repository:
    `$ clone https://github.com/mupq/pqriscv-vexriscv.git`


### How to reproduce the results
First, go to the folder where the thesis project is cloned and compile all levels of Kyber using the following commands:
```sh
$ ./build everything.py -s pqvexriscvsim kyber512
$ ./build everything.py -s pqvexriscvsim kyber768
$ ./build everything.py -s pqvexriscvsim kyber1024
```
All binaries are generated in the `bin` folder and they need to be manually run one by one. In order to reproduce the results, switch to the folder where vexriscv is cloned. We are using the following template command to explain how to run all the binaries:
```sh
sbt "runMain mupq.PQVexRiscvSim --init ../pqriscv/bin/crypto_kem_kyber512_kyber512r1_speed.bin"
```
In this template, `crypto_kem_kyber512_kyber512r1_speed.bin` is the name of the binary that we are executing (generated previously):
Replace **kyber512** with **kyber768** or **kyber1024** to get the results specific to the other security levels of Kyber.
Replace **r1** with **r2** to get the results from round 2, and with **ref1** or **ref2** to get the results from the reference code of round 1 and 2. Each test has to be run individually. The current command produces the speed results indicated by the last word in it: **speed**.
To compute how much Keccak takes in each block (key generation, encapsulation and decapsulation), replace the **speed** word with **hashing**.
To get the **code size** of Kyber round 1 and 2, go to the the thesis folder and then run:
```sh
$ cd crypto kem/code size/round1
$ sh code size.sh
$ cd ../round2
$ sh code size.sh
```
