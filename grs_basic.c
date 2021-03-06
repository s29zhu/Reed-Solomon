/*Elimination Method
*Author:Shasha Zhu
*Jan 2013
*This file contains basic functions of GRS. Including computing the column
*multipliers of parity chechk matrix, computing parity check matrix, given
*the information word, get corresponding codeword etc.
*/

#include "stdio.h"
#include "time.h"
#include "stdlib.h"
#include "galois.h"
#include "grs.h"
#define ERROR_NUM 1

//compute multipliers withe elimination method
/*

in this case, the linear equations are given as a format of

m*n matrix          n*1 matrix  
|a1,a2,a3,...,an|   |x1  |      
|b1,b2,b3,...,bn|   |x2  |
    ...     *   |... |  = 0
|...............|   |xn  |
*/
int * ComputeMultipliers(int code_locators[M][N], 
            int *multipliers){
    //i and k are row index, j is column index
    int i = 0, j = 0, k = 0;
    int add_vector[N] = {0};
    int elimination = 0;
    //printf("****%d\n", galois_single_divide(11, 6, 4)); 
    //eliminate the first M - 1 elements
    for(i = 0; i < M - 1 ; i++){
        elimination = code_locators[i][i];
        for(j = i; j < N; j++){
            //code_locators[i][j] = 0;
            code_locators[i][j] = galois_single_divide(code_locators[i][j],
                                                       elimination,
                                                        w); 
            /*the will be eliminated element is replaced by the add_vector
            multiply by other elements*/
            add_vector[j] = code_locators[i][j];
        }
        /*for every row of the matrix, replace the will-be-eliminated element with
        *the add vector. code_locators[k][i] is the will-be-eliminated element.
        */
        for(k = i + 1; k < M; k++){
        //  code_locators[k][i] = 0;
            elimination = code_locators[k][i];
            for(j = i; j < N; j++){
                code_locators[k][j] ^= galois_single_multiply(add_vector[j],
                                                                elimination,
                                                                w);
            }           
        }
    }

    elimination = code_locators[i][i];
    for(j = i; j < N; j++){
        //code_locators[i][j] = 0;
        code_locators[i][j] = galois_single_divide(code_locators[i][j],
                                                   elimination,
                                                    w);
    } 
    
    /* now there are M equations and there are N multipliers. we set the multiplier[N-1] = 1
    *and hence get the value of other N - 1 i.e. M multipliers.
    */      
    multipliers[N - 1] = 1;
    for(i = N - 2; i >= 0; i--){
        multipliers[i] = 0;
        for(j = i + 1; j < N; j++)
            multipliers[i] ^= galois_single_multiply(code_locators[i][j], multipliers[j], w);
    }
    return multipliers;

}

/* Get the codeword by evaluating the information polynomial
*information polynomial's coefficients are stored in array information from low
degree to high degre
*/
int* ComputeCodeword(int * information, int info_size, int *codeword, int codeword_size){
    int i = 0, j = 0, k = 0;
    int temp_power = 1;
    for(i = 1; i <= codeword_size; i++){
        codeword[i - 1] = 0;
        for(j = 0; j < info_size; j++){
            temp_power = 1;
            for(k = 0; k < j; k++) 
                temp_power = galois_single_multiply(temp_power, i, w);
            codeword[i - 1] ^= galois_single_multiply(information[j], temp_power, w);
        }
    }
    
    return codeword;
}

/* Here we assume the column multiplers are all ones. Hence the code locators matrix is 
*exactly the generator matriix
*/
void GeneratorMatrix(int row_size, int column_size, int (*generator)[column_size]){
    int i = 0, j = 0, k = 0;
    int temp = 1;

    for(i = 0; i < row_size; i++){
        for(j = 0; j < column_size; j++){
            // j + 1 will range from 1 to 8, those are our code locators
            temp = 1;
            // this for loop computes the i-power of every code locators
            for(k = 0; k < i; k++)
                temp = galois_single_multiply(temp, j + 1, w);
            generator[i][j] = temp;
        }
    }
    
}

/*multiply the code locators matrix and column multiplier matrix.
*column multiplier matrix is a diagnal matrix, so here use 1-dimension vector instead 
* for the convenience of programming
*/
void MatrixMultiply(int row_size, int column_size, int (*matrix)[column_size], int *column_multipliers){
    int i = 0, j = 0, k = 0;
    int temp = 0;
    for(i = 0; i < row_size; i++){
        for(j = 0; j < column_size; j++){
            matrix[i][j] = galois_single_multiply(matrix[i][j], column_multipliers[j], w);  
        }
    }
}

/* size of parity check matrix is (N - K)*N, the codeword size is N*1
*so the sydrome is a (N - K)*1 vector
*/

void ComputeParityCheckMatrix(int parity_check_matrix[N - K][N]){
    int multipliers[N] = {0};
    int pre_parity_check_matrix[M][N] = {0};
    int i = 0, j = 0, k = 0;
    

    GeneratorMatrix(M, N, pre_parity_check_matrix);
//copy pre_parity_check_matirx to formal parity_check_matrix because later when
//we call ComputeMultipliers() function, the pre_parity_check_matrix will be
//changed
    for(i = 0; i < N - K; i++){
        for(j = 0; j < N; j++){
            parity_check_matrix[i][j] = pre_parity_check_matrix[i][j];
        }
    }
/*
*Use pre_parity_check_matrix to compute the column multipliers of the parity check matrix
*with elimination method
*/
//compute the colomn multipliers of parity check matrix
// M*N  **   N*1
    ComputeMultipliers(pre_parity_check_matrix, multipliers);
//compute the parity check matrix
    MatrixMultiply(N - K, N, parity_check_matrix, multipliers);

}

void ComputeSyndrome(int row_size, int column_size, 
                        int (*parity)[column_size], 
                        int *codeword,
                        int *syndrome){
    int i = 0, j = 0 , k = 0;
    for(i = 0; i < row_size; i++){
        syndrome[i] = 0;
        for(j = 0;j < column_size; j++){
            syndrome[i] ^= galois_single_multiply(parity[i][j], codeword[j], w);
        }    
    }
}

int main(void){
    
    int information[K] = {4, 9, 13, 7, 14};
    int codeword[N] = {0};
    int multipliers[N] = {0};
    int generator[K][N] = {0};
    int pre_parity_check_matrix[M][N] = {0};
    int parity_check_matrix[N - K][N] = {0};
    int syndrome[N] = {0};
    int i = 0, j = 0, k = 0;
/*
*Use pre_parity_check_matrix to compute the column multipliers of the parity check matrix
*with elimination method
*/
    GeneratorMatrix(M, N, pre_parity_check_matrix);
//copy pre_parity_check_matirx to formal parity_check_matrix
    for(i = 0; i < N - K; i++){
        for(j = 0; j < N; j++){
            parity_check_matrix[i][j] = pre_parity_check_matrix[i][j];
        }
    }
/*
*Use pre_parity_check_matrix to compute the column multipliers of the parity check matrix
*with elimination method
*/
//compute the colomn multipliers of parity check matrix
    ComputeMultipliers(pre_parity_check_matrix, multipliers);
//compute the parity check matrix
    MatrixMultiply(N - K, N, parity_check_matrix, multipliers);
/*compute the generator matrix, since the colomn multipliers of generator matrix are ones,
*then the we eliminate the step to multiply the generator matrix by colomn multipliers.
*/
    GeneratorMatrix(K, N, generator);   
//compute the codeword  
    ComputeCodeword(information, K, codeword, N);
//Assume that there is an error in the codeword, and suppose its on the second
//position
    codeword[2] = 0;
    ComputeSyndrome(N - K, N, parity_check_matrix, codeword, syndrome);
    // here we can see that the syndrom is not all zero vector any more
     unsigned int error_locator_poly[(N - K)/2 + 1] = {0};
     unsigned int error_locator_poly_derivative[(N - K)/2 + 1] = {0};
     unsigned int error_evaluator_poly[N] = {0};
     unsigned int error_magnitudes[ERROR_NUM] = {0};
     unsigned int locators[ERROR_NUM] = {0}; 
     int error_location = 0;
     int *galois_ilog_table = NULL;

    galois_ilog_table = galois_get_ilog_table(w);
     locators[0] = galois_ilog_table[2];
     get_error_locate_poly(locators, ERROR_NUM, error_locator_poly);
     get_evaluator_poly(syndrome, error_locator_poly, error_evaluator_poly);
     get_formal_derivation(error_locator_poly, error_locator_poly_derivative, (N - K)/2 + 1);
     evaluate_error(error_magnitudes, 
                        error_locator_poly_derivative, 
                        error_evaluator_poly, 
                        locators,
                        ERROR_NUM);
     //add the error magnitudes to received codeword
     for(i = 0; i < ERROR_NUM; i++){
         error_location = galois_log(locators[i], w);
         codeword[error_location] ^= error_magnitudes[i];
     }

    return 1;
}
