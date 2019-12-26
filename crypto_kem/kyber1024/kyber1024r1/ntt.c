#include "inttypes.h"
#include "ntt.h"
#include "params.h"
#include "reduce.h"
//mine
#include <stdio.h> 


extern const uint16_t omegas_inv_bitrev_montgomery[];
extern const uint16_t psis_inv_montgomery[];
extern const uint16_t zetas[];

extern void ntt_merged(uint16_t *p,const uint16_t *zetas, uint32_t KQ, uint32_t KQ4,uint32_t qinv,uint16_t *new_zeta_30);

extern void invntt_merged(uint16_t *p,const uint16_t *omegas, uint32_t KQ, uint32_t KQ4,uint32_t qinv,const uint16_t* psis_inv_mont, uint16_t *new_omega_30);

//ROUND 1 zetas for asm fwNTT level 3-0 
// generated from python_scripts/gen_zetas.py
uint16_t new_zeta_30[240] =
{3583, 606, 729, 3465, 6792, 1538, 4664, 7613, 7295, 2007, 323, 5112, 3716, 2289, 6442, 7010, 5383, 962, 2023, 7643, 3660, 7673, 6965, 2713, 7126, 3401, 963, 6596, 607, 5027, 6414, 3240, 7548, 1694, 6905, 3995, 3475, 7078, 4484, 5937, 944, 2860, 2680, 5049, 1777, 263, 5129, 7653, 5939, 1859, 6910, 4434, 5850, 3387, 6487, 6777, 4812, 4724, 7077, 186, 1285, 5929, 4965, 1019, 1492, 7087, 4761, 6848, 6793, 3463, 5877, 1174, 7116, 3077, 5945, 291, 2461, 641, 657, 4859, 5798, 2640, 6591, 590, 6643, 1337, 6036, 3991, 1675, 2053, 7143, 1584, 2666, 1693, 2607, 2782, 5400, 6055, 1162, 1679, 3883, 4311, 2106, 6163, 4486, 7338, 1142, 157, 6466, 1010, 957, 3851, 6374, 5006, 4576, 4288, 5180, 4102, 282, 6119, 1581, 7407, 5222, 2121, 6392, 7319, 3367, 7443, 6330, 3184, 4971, 2530, 5325, 4171, 7185, 5134, 5602, 5142, 3659, 3375, 6430, 7583, 5175, 5655, 1898, 382, 7211, 43, 5965, 6073, 5184, 6140, 5485, 1549, 5856, 4773, 6084, 1730, 332, 1577, 3304, 2329, 1699, 6150, 2379, 5932, 4931, 1559, 5544, 1650, 3997, 4390, 5113, 333, 3502, 4517, 1480, 1172, 5567, 651, 4042, 2085, 5284, 6722, 2915, 4245, 2635, 925, 4573, 599, 1367, 4109, 1863, 6929, 1605, 5775, 2056, 3538, 6128, 7676, 5737, 1616, 3866, 2065, 4048, 839, 5764, 2447, 2022, 3345, 2468, 7269, 3535, 3457, 3132, 7196, 4702, 1990, 4067, 2036, 2069, 3567, 7371, 2368, 339, 3, 7190, 1957, 6239, 851, 2122, 3009, 6947, 2159, 654, 7327, 2768, 6676, 987, 2214};

uint16_t new_omega_30[240] =
{990, 254, 862, 5047, 6586, 5538, 4400, 7103, 990, 254, 862, 5047, 990, 254, 990, 2025, 6804, 3858, 1595, 2299, 4345, 1319, 7197, 6586, 5538, 4400, 7103, 862, 5047, 254, 7678, 5213, 1906, 3639, 1749, 2497, 2547, 6100, 2025, 6804, 3858, 1595, 6586, 5538, 862, 343, 538, 7390, 6396, 7418, 1267, 671, 4098, 2299, 4345, 1319, 7197, 4400, 7103, 5047, 5724, 491, 4146, 412, 4143, 5625, 2397, 5596, 7678, 5213, 1906, 3639, 2025, 6804, 6586, 6122, 2750, 2196, 1541, 2539, 2079, 2459, 274, 1749, 2497, 2547, 6100, 3858, 1595, 5538, 7524, 6539, 5015, 6097, 7040, 5220, 2716, 1752, 343, 538, 7390, 6396, 2299, 4345, 4400, 28, 2552, 133, 4441, 6719, 2298, 6952, 7075, 7418, 1267, 671, 4098, 1319, 7197, 7103, 4672, 5559, 6830, 1442, 2979, 485, 4549, 4224, 5724, 491, 4146, 412, 7678, 5213, 2025, 6065, 1944, 5, 1553, 5046, 3436, 4766, 959, 4143, 5625, 2397, 5596, 1906, 3639, 6804, 3291, 3684, 6031, 2137, 1597, 2908, 1825, 6132, 6122, 2750, 2196, 1541, 1749, 2497, 3858, 98, 1251, 4306, 4022, 4314, 362, 1289, 5560, 2539, 2079, 2459, 274, 2547, 6100, 1595, 3830, 6724, 6671, 1215, 2281, 4899, 5074, 5988, 7524, 6539, 5015, 6097, 343, 538, 2299, 5041, 1883, 2822, 7024, 2920, 594, 6189, 6662, 7040, 5220, 2716, 1752, 7390, 6396, 4345, 3247, 771, 5822, 1742, 4206, 3686, 776, 5987, 28, 2552, 133, 4441, 7418, 1267, 1319, 8, 4021, 38, 5658, 3017, 6143, 889, 4216, 6719, 2298, 6952, 7075, 671, 4098, 7197};


/*************************************************
* Name:        ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial (vector of 256 coefficients) in place;
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint16_t *p: pointer to in/output polynomial
**************************************************/
void ntt(uint16_t *p)
{
  /* ASM function. Merged level 7-4 and 3-0 */
 ntt_merged(p,zetas, KYBER_Q, KYBER_Q<<2, 7679U,new_zeta_30); //p, z, Q, 4*Q, qinv =7679  


 /*int level, start, j, k;
  uint16_t zeta, t;

  k = 1;
  for(level = 7; level >= 0; level--)
  {
    for(start = 0; start < KYBER_N; start = j + (1<<level))
    {
      zeta = zetas[k++];
      for(j = start; j < start + (1<<level); ++j)
      {
        t = montgomery_reduce((uint32_t)zeta * p[j + (1<<level)]);

        p[j + (1<<level)] = barrett_reduce(p[j] + 4*KYBER_Q - t);
           
        if(level & 1) // odd level 
          p[j] = p[j] + t; // Omit reduction (be lazy) 
        else         
          p[j] = barrett_reduce(p[j] + t);
          
      }
    }
  }*/
  
}

  
/*************************************************
* Name:        invntt
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT) of
*              a polynomial (vector of 256 coefficients) in place;
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint16_t *a: pointer to in/output polynomial
**************************************************/
void invntt(uint16_t * a)
{
  /* ASM function. Merged level 0-3 and 4-7 */
  invntt_merged(a,omegas_inv_bitrev_montgomery, KYBER_Q, KYBER_Q<<2, 7679U,psis_inv_montgomery,new_omega_30); //qinv =7679 


 /* int start, j, jTwiddle, level;
  uint16_t temp, W;
  uint32_t t;
 
  for(level=0;level<8;level++)
  {
   
    for(start = 0; start < (1<<level);start++)
    {
      jTwiddle = 0;

      for(j=start;j<KYBER_N-1;j+=2*(1<<level))
      {
        W = omegas_inv_bitrev_montgomery[jTwiddle++];              
        temp = a[j];

        if(level & 1) // odd level 
          a[j] = barrett_reduce((temp + a[j + (1<<level)]));           
        else
          a[j] = (temp + a[j + (1<<level)]); // Omit reduction (be lazy) 
          
        t = (W * ((uint32_t)temp + 4*KYBER_Q - a[j + (1<<level)]));

        a[j + (1<<level)] = montgomery_reduce(t);        
      }
    }
  }

  for(j = 0; j < KYBER_N; j++)
    a[j] = montgomery_reduce((a[j] * psis_inv_montgomery[j]));*/
}


