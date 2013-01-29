/*
*Author:Shasha Zhu
*Jan 2013
*
*In the is file, consumer get secret shares in the form of e(g,g)^(s1*s), si*is a secret share. In this case, i ranges from 0 to 7 (since we add two dum
*my nodes to the tree, the number of shares now is 8 instead of 6). Si is ai*lso q(i). q(x) is root polynomial. See grs_encode_check.c for details.
*/

#include <pbc.h>
#include <pbc_test.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include "grs.h"

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

void consumerShares(signed long int *codeword){
    pairing_t pairing;
    element_t g, r, a, e_g_g, share;
    char *argv = "./param/a.param";
    char s[16384];
    signed long int temp_share;
    FILE *fp = stdin;

    fp = fopen(argv, "r");
    if (!fp) 
        pbc_die("error opening %s\n", argv);
    size_t count = fread(s, 1, 16384, fp);
    if(!count) 
        pbc_die("read parameter failure\n");
    fclose(fp);
    if(pairing_init_set_buf(pairing, s, count)) 
        pbc_die("pairing init failed\n");
    if(!pairing_is_symmetric(pairing)) pbc_die("pairing is not symmetric\n");
    
    element_init_G1(g, pairing);
    element_init_Zr(r, pairing);
    element_init_Zr(a, pairing);
    element_init_Zr(share, pairing);
    element_init_GT(e_g_g, pairing);
    
    //find the generator of the group
    element_set(g, ((curve_data_ptr)((a_pairing_data_ptr)
    pairing->data)->Eq->data)->gen);
    element_random(r);
    element_random(a);
    //compute e(g, g)
    element_pairing(e_g_g, g, g);
    //compute e(g, g)^r
    element_pow_zn(e_g_g, e_g_g, r);
    //compute e(g,g)^ra
    element_pow_zn(e_g_g, e_g_g, a);
    temp_share = codeword[0];
    //transfer signed long int type ecret shares to an element_t type before we do the power of
    //e_g_g
    element_set_si(share, temp_share);
    element_pow_zn(e_g_g, e_g_g, share);
    
}
