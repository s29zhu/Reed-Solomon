#include <pbc.h>
#include <pbc_test.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <string.h>

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
    element_t g, alpha, e_g_g, secret_key, pub_key;
    FILE *fp;
    char s[16384];

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
    element_init_Zr(alpha, pairing);
    element_init_G1(secret_key, pairing);
    element_init_GT(pub_key, pairing);
    

    element_set(g, ((curve_data_ptr)((a_pairing_data_ptr)
    pairing->data)->Eq->data)->gen);
    element_random(alpha);
    element_pow_zn(secret_key, g, alpha);
    element_pairing(pub_key, g, secret_key);
    fp = fopen("owner_secret_key.txt", "w+");
    if(!fp)
        pbc_die("error opening owner secret key file.");
    else{
        fprintf(fp, "secret key:");
        element_out_str(fp, 10, secret_key);
        fclose(fp);
    }
    fp = fopen("../public/owner_public_key.txt", "w+");
    if(!fp)
        pbc_die("error opening error public key file.");
    else{
        fprintf(fp, "public key:");
        element_out_str(fp, 10, pub_key);
        fclose(fp);
    }
    element_clear(g);
    element_clear(alpha);
    element_clear(secret_key);
    element_clear(pub_key);
    return 1;
}
