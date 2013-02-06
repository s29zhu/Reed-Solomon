#include <stdio.h>
#include <stdlib.h>
#include "galois.h"
#include "polynomial.h"
#include <time.h>

//This function is used to compute (1+(A^i)x)(1+(A^k)x)...., where i and k are the location of where went wrong
//In this case, we put the secret share as at the fixed location and thus we already know about the location of errors
void get_error_locate_poly(unsigned int *locators, int length_of_locators, unsigned int *error_locator_poly){
    int i = 0, j = 0;   
    unsigned int pre_error_locator_poly[N - K + 1] = {0};
    
    error_locator_poly[N - K] = 1;
    for(i = 0; i < length_of_locators; i++){
        for(j = 0; j < N - K + 1; j++){
            pre_error_locator_poly[j] = error_locator_poly[j];
        }
        for(j = N - K; j >= 0 ; j--){
            error_locator_poly[j] = galois_single_multiply(error_locator_poly[j], locators[i], w);                  
        }
        for(j = 0; j < N - K + 1 ; j++){
            error_locator_poly[j] = error_locator_poly[j + 1];                  
        }
        error_locator_poly[N - K] = 0;
        for(j = 0; j < N - K + 1; j++){
            error_locator_poly[j] ^= pre_error_locator_poly[j];
        }
    }
}

// compute the derivative of given polynomial poly
void formal_derivative(unsigned int *poly, unsigned int *derivative, int poly_length){
    int i = 0;
    for(i = poly_length - 1; i >= 0; i--){
        if((N - K - i - 1) % 2 == 0)
            derivative[i] = 0;  
        else
            derivative[i] = poly[i];
    }
    for(i = poly_length; i > 0; i--){
        poly[i] = poly[i - 1];
    }
    poly[0] = 0;
}

void get_evaluator_poly(unsigned int *syndrome, unsigned int *evaluator_poly){
    
}

void evaluate_error(unsigned int ){

}

//count the non-zero elements of array
void get_length_Of_Array(unsigned int *array, unsigned int *length){
    unsigned int ele = *array;
    *length = (1 << w) - 1;
    while(!ele){
        *length -= 1;
        array += 1;
        ele = *array;   
    }
}

//divident polynomial divided by divisor. The reminder is stored in divident. 
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
/*  printf("The remainder is:\n");
    for(i; i < (1 << w) - 1; i++)
        printf("%d  ", divident[i]);
    printf("\n");*/
}

int main(void){
    unsigned int poly[4] = {5, 7, 3, 1};//A^1, A^4
    unsigned int derivative[7] = {0};
    unsigned int error_locator_poly[N - K + 1] = {0};
    int i = 0;

//  convolute_poly(locators, 3, error_locator_poly);
    formal_derivative(poly, derivative, 4);
    for(i = 0; i < N - K + 1; i++){
        printf("%d  ", derivative[i]);
    }

    printf("\n");
    return 1;
}
