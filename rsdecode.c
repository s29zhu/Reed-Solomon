/*rsdecode.c
 *Shasha Zhu
 *Oct, 2012

The decode algorithm is an implementation of page219 of Error Control Systems 
For Digital Communication and Storage. 

In general case the Berlekamp_Massey_decode() function is used to compute the 
error locator polynomial, then we need to factorize this polynomial to get the 
error locations. In this secret sharing case, we assume that we already know the
error locations and then we can build error locator polynomial use
get_error_locate_poly() function from the polynomial.c file.
*/
#include "stdio.h"
#include "stdlib.h"
#include "galois.h"
#include "polynomial.h"
#define ERROR_NUM 1

static unsigned int parity_check_matrix[N - K][N] = {0};

void get_parity_check_matrix(){
	int *galois_ilog_table = NULL;
	int i = 0, j = 0;
	galois_ilog_table = galois_get_ilog_table(w);
	for(j = 0; j < N - K; j++){
		if(j == 0)
			for(i = 0; i < N; i++)
				parity_check_matrix[j][i] = galois_ilog_table[i];
		else
			for(i = 0; i < N; i++)
				parity_check_matrix[j][i] = galois_single_multiply(parity_check_matrix[j - 1][i], galois_ilog_table[i], w);	
	}
/*	printf("The parity check matrix is:\n");
	for(i = 0; i < N - K; i++){
		for(j = 0; j < N; j++){
			printf("%d	", parity_check_matrix[i][j]);
		}
		printf("\n");
	}*/
}

//when compute the syndrome, codeword should be in format of c0c1c2..., other
//than the polynomial reverse format
//Syndrome = [s1, s2, ..., sr]
unsigned int *compute_syndrome(unsigned int *syndrome, unsigned int *received_codeword){
	unsigned int sum = 0;
	int i = 0, j = 0;
	
	get_parity_check_matrix();
	for(i = 0; i< N - K; i++){
		sum = 0;
		for(j = 0; j < N; j++){
			sum ^= galois_single_multiply(parity_check_matrix[i][j], received_codeword[j], w);
		}
		syndrome[i] = sum;
	}
	printf("The syndrome is:\n");
	for(i = 0; i < N - K; i++)
		printf("%d	", syndrome[i]);
	printf("\n");
	return syndrome;
}

//
/*void Berlekamp_Massey_decode(unsigned int *syndrome, unsigned int *locators){
	int k = 0, L = 0;
	unsigned int error_locator_poly[N - K + 1] = {0}; 
	unsigned int temp_poly[N - K + 1] = {0};
	unsigned int discrepancy = 0;
	unsigned int multi_temp = 0
	
	error_locator_poly[N - K] = 1;
	temp_poly[N - K - 1] = 1;
	
	for(k = 1; k < N - K; k++){
		multi_temp = 0;
		for(i = 1; i <= L; i++)
			multi_temp ^= galois_single_multiply(error_locator_poly[N - K - i], syndrome[N - 1 - (k - i)], w);
		discrepancy = syndrome[N - k - 1]^multi_temp;
		if(discrepancy == 0){			
			break;
		}
		for(i = 0; i < N - K + 1; i++){
			error_locator_poly[i] ^= galois_single_multiply(discrepancy, temp_poly[i], w);
		}
		if(2*L < k){
			L = k - L;
			for(i = 0; i < N - K + 1; i++){
				temp_poly[i] = galois_single_divide(error_locator_poly[i], discrepancy, w);
			}
		}
		for(i = 0; i < N - K; i++){
			temp_poly[i] = temp_poly[i + 1];
		}
		temp_poly[N - K] = 0;
	}
	// add factorization code here
	factorize_poly(locators, error_locator_poly);///????hard problem
}*/

//
void evaluate_error(unsigned int *error_magnitudes, 
			unsigned int *error_locator_poly_derivative, 
			unsigned int *error_evaluator_poly, 
			unsigned int *locators, 
			int length_of_locators){
	int i = 0, j = 0, k = 0;
	unsigned int poly_evaluation = 0;
	unsigned int power = 1;
	unsigned int inverse_locator = 0 ;
	for(i = 0; i < length_of_locators; i++){
		inverse_locator = galois_inverse(locators[i], w);
		//nominator computing
		poly_evaluation = 0;		
		for(j = N - 1; j >= K - 1; j--){//?
			power = 1;
			for(k = 0; k < N - 1 - j; k++)
				power = galois_single_multiply(power, inverse_locator, w);
			poly_evaluation ^= galois_single_multiply(power, error_evaluator_poly[j], w); 
		//	printf("test %d	\n", poly_evaluation);
		}		
		error_magnitudes[i] = galois_single_multiply(poly_evaluation, locators[i], w);
		//denominator computing
		poly_evaluation = 0;
		for(j = (N - K)/2; j >= 0; j--){
			power = 1;
			for(k = 0; k < (N - K)/2 - j; k++)
				power = galois_single_multiply(power, inverse_locator, w);
			poly_evaluation ^= galois_single_multiply(power, error_locator_poly_derivative[j], w); 
		}
		// error magnitudes  computing		
		error_magnitudes[i] = galois_single_divide(error_magnitudes[i], poly_evaluation, w);
	}
	printf("The error magnitudes are:\n");
	for(i = 0; i < length_of_locators; i ++){
		printf("%d	", error_magnitudes[i]);
	}
	printf("\n");
}

void rsdecode(unsigned int *locators, unsigned int *received_codeword){
	unsigned int syndrome[N - K] = {0};
	unsigned int error_locator_poly[(N - K)/2 + 1] = {0};
	unsigned int error_locator_poly_derivative[(N - K)/2 + 1] = {0};
	unsigned int error_evaluator_poly[N] = {0};
	unsigned int error_magnitudes[ERROR_NUM] = {0};
	int i = 0;
	int error_location = 0;

	compute_syndrome(syndrome, received_codeword);
	get_error_locate_poly(locators, ERROR_NUM, error_locator_poly);
	get_evaluator_poly(syndrome, error_locator_poly, error_evaluator_poly);
	get_formal_derivation(error_locator_poly, error_locator_poly_derivative, (N - K)/2 + 1);
	evaluate_error(error_magnitudes, error_locator_poly_derivative, error_evaluator_poly, locators,ERROR_NUM);
	//add the error magnitudes to received codeword
	for(i = 0; i < ERROR_NUM; i++){
		error_location = galois_log(locators[i], w);
		received_codeword[error_location] ^= error_magnitudes[i];		
	}
}
/*
int main(void){
	int i = 0; 
	//error_locator_poly is in a format of (1 + (A^i)x)*...*(1 + (A^j)x)
	//locators[i] represents the  A^i, ..., A^j in the error locator polynomial
	unsigned int locators[ERROR_NUM] = {6};
	unsigned int received_codeword[N]={2, 7, 3, 1, 0, 6, 4};
	rsdecode(locators, received_codeword);
	printf("After correction, the codeword is:\n");
	for(i = 0; i < N; i++){
		printf("%d	", received_codeword[i]);
	}
	printf("\n");
	return 1;
}*/
