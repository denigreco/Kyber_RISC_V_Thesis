#include <stdint.h>
#include <string.h>

#include "address.h"
#include "params.h"
#include "thash.h"

#include "fips202.h"

/**
 * Takes an array of inblocks concatenated arrays of SPX_N bytes.
 */
static void PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash(
    unsigned char *out, unsigned char *buf,
    const unsigned char *in, unsigned int inblocks,
    const unsigned char *pub_seed, uint32_t addr[8]) {

    unsigned char *bitmask = buf + SPX_N + SPX_ADDR_BYTES;
    unsigned int i;

    memcpy(buf, pub_seed, SPX_N);
    PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_addr_to_bytes(buf + SPX_N, addr);

    shake256(bitmask, inblocks * SPX_N, buf, SPX_N + SPX_ADDR_BYTES);

    for (i = 0; i < inblocks * SPX_N; i++) {
        buf[SPX_N + SPX_ADDR_BYTES + i] = in[i] ^ bitmask[i];
    }

    shake256(out, SPX_N, buf, SPX_N + SPX_ADDR_BYTES + inblocks * SPX_N);
}

/* The wrappers below ensure that we use fixed-size buffers on the stack */

void PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash_1(
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]) {

    unsigned char buf[SPX_N + SPX_ADDR_BYTES + 1 * SPX_N];
    PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash(
        out, buf, in, 1, pub_seed, addr);
}

void PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash_2(
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]) {

    unsigned char buf[SPX_N + SPX_ADDR_BYTES + 2 * SPX_N];
    PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash(
        out, buf, in, 2, pub_seed, addr);
}

void PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash_WOTS_LEN(
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]) {

    unsigned char buf[SPX_N + SPX_ADDR_BYTES + SPX_WOTS_LEN * SPX_N];
    PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash(
        out, buf, in, SPX_WOTS_LEN, pub_seed, addr);
}

void PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash_FORS_TREES(
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]) {

    unsigned char buf[SPX_N + SPX_ADDR_BYTES + SPX_FORS_TREES * SPX_N];
    PQCLEAN_SPHINCSSHAKE256128FROBUST_CLEAN_thash(
        out, buf, in, SPX_FORS_TREES, pub_seed, addr);
}
