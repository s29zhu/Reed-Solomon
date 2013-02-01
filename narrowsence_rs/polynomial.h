//RS w-bit per symbol
#define w 3 
//RS codeword length
#define N (i<<w - 1)
#define M N - 1
//#define Num  60
//RS information word length
#define K 3

void factorization(unsigned int *poly);
void lengthOfArray(unsigned int *array, unsigned int *length);
void remainderComputing(unsigned int *devident, unsigned int *devisor);
void get_evaluator_poly(unsigned int *syndrome, unsigned int *error_locator_poly, unsigned int *error_evaluator_poly);
void get_formal_derivation(unsigned int *poly, unsigned int *derivative, int poly_length);
void get_error_locate_poly(unsigned int *locators, int length_of_locators, unsigned int *error_locator_poly);
