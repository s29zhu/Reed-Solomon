/*polynomial.c
 *Shasha Zhu
 *Oct, 2012

This file contains several functions that related to RS encoding 
and decoding polynomial computing.
note: We use array to represent the polynomial. And the format 
is different from codeword. For example, poly[5] = {3,4,5,7,1} means 
1 + 7x + 5x^2 + 4x^3 + 3x^4.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.	
*/

#include <stdio.h>
#include <stdlib.h>
#include "galois.h"
#include "polynomial.h"
#include <time.h>

//This function is used to compute (1+(A^i)x)(1+(A^k)x)...., where i and k are the location of where went wrong
//In this case, we put the secret share as at the fixed location and thus we already know about the location of errors
void get_error_locate_poly(unsigned int *locators, int length_of_locators, unsigned int *error_locator_poly){
	int i = 0, j = 0;	
	unsigned int pre_error_locator_poly[(N - K)/2 + 1] = {0};
	
	error_locator_poly[(N - K)/2] = 1;
	for(i = 0; i < length_of_locators; i++){
		//save the re-multiplication polynomial
		for(j = 0; j < (N - K)/2 + 1; j++){
			pre_error_locator_poly[j] = error_locator_poly[j];
		}
		//multiply the locators one by one
		for(j = (N - K)/2; j >= 0 ; j--){
			error_locator_poly[j] = galois_single_multiply(error_locator_poly[j], locators[i], w); 					
		}
		//move the multiplication result to one higher degree
		for(j = 0; j < (N - K)/2 + 1 ; j++){
			error_locator_poly[j] = error_locator_poly[j + 1]; 					
		}
		error_locator_poly[(N - K)/2] = 0;
		//added to the revious error_locator_poly
		for(j = 0; j < (N - K)/2 + 1; j++){
			error_locator_poly[j] ^= pre_error_locator_poly[j];
		}
	}
	printf("error_locator_poly:\n");
	for(i = 0; i < (N - K)/2 + 1; i++){
		printf("%d	", error_locator_poly[i]);
	}
	printf("\n");
}

// get formal derivation of the given poly
void get_formal_derivation(unsigned int *poly, unsigned int *derivative, int poly_length){
	int i = 0;
	for(i = poly_length - 1; i >= 0; i--){
		if((poly_length - i - 1) % 2 == 0)
			derivative[i] = 0;	
		else
			derivative[i] = poly[i];
	}
	for(i = poly_length - 1; i > 0; i--){
		derivative[i] = derivative[i - 1];
	}
	derivative[0] = 0;
	printf("The derivative is: \n");
	for(i = 0; i < poly_length; i++){
		printf("%d	", derivative[i]);
	}
	printf("\n");
}

//size of syndrome is r, size of error_locator_poly is (N - K) / 2 + 1, size of error_evaluator_poly is N
//evaluator_poly = (1 + syndrome_poly)*(error_locator_poly)
//This function only works with berley-kalm algorithm for narrow sense RS

void get_evaluator_poly(unsigned int *syndrome, unsigned int *error_locator_poly, unsigned int *error_evaluator_poly){
	int i = 0 , j = 0, k = 0, l = 0;
	unsigned int syndrome_plus_one[N - K + 1] = {0};
	//compute 1 + s(x), s(x) is the syndrome polynomial
	for(i = 1 ; i < N - K + 1; i++){
		syndrome_plus_one[i] = syndrome[i - 1];
	}
	syndrome_plus_one[0] = 1;
	//compute (1 + s(x))*error_locator_poly
	for(i = 0; i <= N - K; i++){ //s1, s2, ..., sr multiply error_locator_poly one after each other
		for(j = (N - K) / 2; j >= 0; j--){
			error_evaluator_poly[N - k - 1 - (N - K) / 2 + j] ^= galois_single_multiply(error_locator_poly[j], syndrome_plus_one[i], w);
			
		}
		k++;
	}  
	//mod x^(2t+1)
	for(k = K - 2; k >= 0; k --)
		error_evaluator_poly[k] = 0;
	printf("The error_evaluator_poly is:\n");
	for(i = 0; i < N; i++)
		printf("%d	", error_evaluator_poly[i]);
	printf("\n");    				
}

//get GRS error evaluate polynomial = segma(j belongs to J)multi(1-alpa(subj)x)
//J is a set of error locations
void grs_get_evaluator_poly(){


}

void get_length_Of_Array(unsigned int *array, unsigned int *length){
	unsigned int ele = *array;
	*length = (1 << w) - 1;
	while(!ele){
		*length -= 1;
		array += 1;
		ele = *array;	
	}
}

void remainderComputing(unsigned int *divident, unsigned int *divisor){
	unsigned int temp = 0;
	unsigned int length_divident = 0;
	unsigned int length_divisor = 0;
	int i = 0;
	// j2 is used to point out the highest degree position of divisor
 	int j2 = 0; 

	while(0 == divisor[j2]){
		j2 += 1;
	}
	get_length_Of_Array(divident, &length_divident);
	get_length_Of_Array(divisor, &length_divisor);
	while(length_divident >= length_divisor){
		int quotient = 0;
		//j1 is used to point out the highest degree position of divident
		int j = 0, j1 = 0;
		while(0 == divident[j1]){	
			j1 += 1;
		}
		quotient = galois_single_divide(divident[j1], divisor[j2], w);
		for(j = j2; j < (1 << w) - 1; j++){
			//sleep(0.1);
			temp = galois_single_multiply(divisor[j], quotient, w);
			divident[j1] = divident[j1] ^ temp;
			j1 ++;
		}		
		get_length_Of_Array(divident, &length_divident);
		get_length_Of_Array(divisor, &length_divisor);
		//printf("length1 %d length2 %d ", length_divident, length_divisor);
	}
	printf("The remainder is:\n");
	for(i; i < (1 << w) - 1; i++)
		printf("%d	", divident[i]);
	printf("\n");
}

