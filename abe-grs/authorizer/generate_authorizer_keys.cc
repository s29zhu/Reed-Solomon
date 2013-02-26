/*Author:Shasha Zhu
*Feb 2013
*/

#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "pbc.h"

typedef struct {
  field_t Fq, Fq2, Eq;
  int exp2, exp1;
  int sign1;
} *a_pairing_data_ptr;

typedef struct {
  field_ptr field; // The field where the curve is defined.
  element_t a, b;  // The curve is E: Y^2 = X^3 + a X + b.
  // cofac == NULL means we're using the whole group of points.
  // otherwise we're working in the subgroup of order #E / cofac,
  // where #E is the number of points in E.
  mpz_ptr cofac;
  // A generator of E.
  element_t gen_no_cofac;
  // A generator of the subgroup.
  element_t gen;
  // A non-NULL quotient_cmp means we are working with the quotient group of
  // order #E / quotient_cmp, and the points are actually coset
  // representatives. Thus for a comparison, we must multiply by quotient_cmp
  // before comparing.
  mpz_ptr quotient_cmp;
} *curve_data_ptr;

int main(void){
    pairing_t pairing;
    element_t g, h, f, beta, beta_inverse;

    char s[16384];
   //signed long int temp_share;
    FILE *fp = stdin;

    fp = fopen("../public/a.param", "r");
    if (!fp) 
        pbc_die("error opening parameter file");
    size_t count = fread(s, 1, 16384, fp);
    if(!count) 
        pbc_die("read parameter failure\n");
    fclose(fp);
    if(pairing_init_set_buf(pairing, s, count)) 
        pbc_die("pairing init failed\n");
    if(!pairing_is_symmetric(pairing)) pbc_die("pairing is not symmetric\n");
    
    element_init_G1(g, pairing);
    element_init_G1(h, pairing);
    element_init_G1(f, pairing);
    element_init_Zr(beta, pairing);
    element_init_Zr(beta_inverse, pairing);

    //(G1, g, h, f) is the public key of authorizer
    //find the generator of the group
    element_set(g, (*(curve_data_ptr) (*(*(a_pairing_data_ptr) (*pairing).data).Eq).data).gen);
    element_random(beta);
    element_invert(beta_inverse, beta);
    //h = g^beta
    element_pow_zn(h, g, beta);
    //f = g^(1/beta)
    element_pow_zn(f, g, beta_inverse);
    fp = NULL;
    fp = fopen("../public/authorizer_public_keys.txt", "w+");
    if(!fp)
        pbc_die("error creating public key files");
    else{
        fprintf(fp, "g:");
        element_out_str(fp, 10, g);
        fprintf(fp, "\n\nh:");
        element_out_str(fp, 10, h);
        fprintf(fp, "\n\nf:");
        element_out_str(fp, 10, f);
        fclose(fp);
    }
    fp = fopen("./authorizer_secret_key.txt", "w+");
    if(!fp)
        pbc_die("error creating secret key files");
    else{
        fprintf(fp, "beta:");
        element_out_str(fp, 10, beta);
    }
    element_clear(g);
    element_clear(h);
    element_clear(f);
    element_clear(beta);
    element_clear(beta_inverse);
    return 1;
}    
    
