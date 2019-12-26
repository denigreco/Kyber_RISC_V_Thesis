#include "api.h"
#include "hal.h"
#include "poly.h"
#include "polyvec.h"
//#include "reduce.h"
//#include "params.h"


#include <stdio.h>
#include <stdint.h>
#include <string.h>

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x####y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_BYTES           NAMESPACE(CRYPTO_BYTES)
#define MUPQ_CRYPTO_PUBLICKEYBYTES  NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES  NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_CIPHERTEXTBYTES NAMESPACE(CRYPTO_CIPHERTEXTBYTES)
#define MUPQ_CRYPTO_ALGNAME NAMESPACE(CRYPTO_ALGNAME)

#define MUPQ_crypto_kem_keypair NAMESPACE(crypto_kem_keypair)
#define MUPQ_crypto_kem_enc NAMESPACE(crypto_kem_enc)
#define MUPQ_crypto_kem_dec NAMESPACE(crypto_kem_dec)


//extern int32_t mont2(int32_t a, uint32_t qinv, uint32_t q);

static void printcycles(const char *s, unsigned long long c)
{
  char outs[32];
  hal_send_str(s);
  snprintf(outs,sizeof(outs),"%llu\n",c);
  hal_send_str(outs);
}

int main(void)
{
  unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
  unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
  unsigned char ct[MUPQ_CRYPTO_CIPHERTEXTBYTES];
  unsigned long long t0, t1;

  hal_setup(CLOCK_BENCHMARK);

  hal_send_str("==========================");

  // Key-pair generation
  t0 = hal_get_time();
  MUPQ_crypto_kem_keypair(pk, sk);
  t1 = hal_get_time();
  printcycles("keypair cycles:", t1-t0);

  // Encapsulation
  t0 = hal_get_time();
  MUPQ_crypto_kem_enc(ct, key_a, pk);
  t1 = hal_get_time();
  printcycles("encaps cycles: ", t1-t0);

  // Decapsulation
  t0 = hal_get_time();
  MUPQ_crypto_kem_dec(key_b, ct, sk);
  t1 = hal_get_time();
  printcycles("decaps cycles: ", t1-t0);

 // NTT
  int16_t poly[256];
  poly_ntt(poly);
  
  t0 = hal_get_time();
  poly_ntt(poly);
  t1 = hal_get_time();
  printcycles("NTT cycles: ",t1-t0);

  // iNTT
  poly_invntt(poly);

  t0 = hal_get_time();
  poly_invntt(poly);
  t1 = hal_get_time();
  printcycles("iNTT: ",t1-t0);



  if (memcmp(key_a, key_b, MUPQ_CRYPTO_BYTES)) {
    hal_send_str("ERROR KEYS\n");
  }
  else {
    hal_send_str("OK KEYS\n");
  }

  hal_send_str("#");
  while(1);
  return 0;
}
