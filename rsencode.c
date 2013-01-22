/*polynomial.c
 *Shasha Zhu
 *Oct, 2012
	
*/

#include "stdio.h"
#include "stdlib.h"
#include "galois.h"
#include "polynomial.h"

unsigned int message[3] = {2, 7, 3};
unsigned int remain[7] = {0};
unsigned int generator[7] = {0, 0, 1, 4, 5, 1, 5};
unsigned int codeword[7] = {0};
unsigned int length = 0;


unsigned int *get_generator_poly();

unsigned int *get_generator_poly(unsigned int *generator_poly){
	unsigned int locators[50] = {0};
	int * galois_ilog_table = NULL;
	int i = 0, j = 0;

	galois_ilog_table = galois_get_ilog_table(w);
	for(i = 0; i < 50; i++){
		locators[i] = galois_ilog_table[i];
	}
	get_error_locate_poly(locators, 50, generator_poly);
}

//format of codeword is (c0, c1, c2, ..., cN-1)
//if we put the codeword in polynomial representation, then it should be like (cN-1, ..., c2, c1, c0)
unsigned int *rsencode(unsigned int *message, unsigned int *generator, unsigned int *codeword){
 	int i = 0;
	for(i = 0; i < K; i++){
		remain[i] = message[i];	
		codeword[N - i - 1] = message[i];
	}
	remainderComputing(remain, generator);
	for(i = 0; i < N - K ; i++){
		codeword[i] ^= remain[N - 1 - i];
	}
	printf("The codeword is:\n");
	for(i = 0; i < N; i++){
		printf("%d	", codeword[i]);
	}
	printf("\n");
	return codeword;
}

int main(void){
	rsencode(message, generator, codeword);
	return 1;
}
