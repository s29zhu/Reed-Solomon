/*rsencode.h
 *Shasha Zhu
 *Dec 2012

This file contains GF(2^w) RS code generators
Note that the prime polynomial representation is different from the representation in galois.c file.
For example, 1 + x^2 + x^3 in galois.c representation is 001011
*/
unsigned int generator3[7] = {0, 0, 1, 4, 5, 1, 5};//w = 3,N = 7, K = 3, prime_poly = 13
unsigned int generator4[15] = {};//w = 4, N = 15, K = 5, prime_poly = 25
unsigned int generator5[31] = {};//w = 5, N = 31, K = 21, prime_poly = 41
unsigned int generator6[63] = {};//w = 6, N = 63, K = 53, prime_poly = 97
unsigned int generator7[127] = {};//w = 7, N = 127, K = 117, prime_poly = 145
unsigned int generator8[255] = {};//w = 8, N = 255, K = 245, prime_poly = 369
unsigned int generator9[511] = {};//w = 9, N = 511, K = 501, prime_poly = 545
unsigned int generator10[1023] = {};//w = 10, N = 1023, K = 1013, prime_poly = 1153
unsigned int generator11[2047] = {};//w = 11, N = 2047, K = 2037, prime_poly = 2561
unsigned int generator12[4095] = {};//w = 12, N = 4095, K = 4085, prime_poly = 6465
unsigned int generator13[8191] = {};//w = 13, N = 8191, K = 8181, prime_poly = 13825
unsigned int generator14[16383] = {};//w = 14, N = 16383, K = 16373, prime_poly = 24849
unsigned int generator15[32767] = {};//w = 15, N = 32767, K = 32757, prime_poly = 49153
unsigned int generator16[65535] = {};//w = 16, N = 65535, K = 65545, prime_poly = 106513
unsigned int generator17[131071] = {};//w = 17, N = 131071, K = 131061
