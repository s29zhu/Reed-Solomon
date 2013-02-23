#include "stdlib.h"
#include "stdio.h"
#include <pbc.h>
#include <pbc_test.h>
#include <math.h>
#include <gmp.h>
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

//qR is a pointer which points to access policy tree
void creat_acc_tree(element_t *qR, pairing_t *pairing);
void init(element_t *qR, element_t ele_att_hash, void *attr[NUM_SHARE]);


void creat_acc_tree(element_t *qR, pairing_t *pairing){
    int i = 0; 

    for(i = 0; i < NUM_SHARE + 1; i++){
        //qR[0] is the secret s
        element_init_Zr(qR[i], *pairing);
        element_random(qR[i]); 
    }
}

void hash_attributes(element_t *ele_attr_hash, void *attr[NUM_SHARE]){
    int i = 0; 
    for(i = 0 ; i < NUM_SHARE + 1; i++){
        element_from_hash(ele_attr_hash[i], attr[i], 256);
    }
}


void initialize(element_t *qR, 
                element_t *ele_att_hash, 
                pairing_t *pairing,
                element_t *g,
                element_t *Cy,
                element_t *Cy_prime){
   //initialize the attribute string 
    int i = 0;
    char *attr[NUM_SHARE];
    FILE *fp = NULL;
    attr[0] = "http://www.servername.com/pathtofile";
    attr[1] = "Julia Zhu";
    attr[2] = "4";
    attr[3] = "r";
    attr[4] = "2013/01/10/13/20";
    attr[5] = "julia@facebook";
    attr[6] = "printer";
    attr[7] = "check string 1";
    attr[8] = "check string 2";
    element_init_G1(*g, *pairing);
    //initialize the attributes' hash value
    for(i = 0; i < NUM_SHARE; i++){
        element_init_G1(ele_att_hash[i], *pairing);
        //element_from_hash(ele_att_hash[i], (void *)attr[i], 256);
        element_init_Zr(qR[i], *pairing);
        element_init_G1(Cy[i], *pairing);
        element_init_G1(Cy_prime[i], *pairing);
    }
    element_init_Zr(qR[i], *pairing);
}

void finalize(element_t *qR,
                element_t *ele_att_hash,
                element_t *g,
                element_t *Cy,
                element_t *Cy_prime){
    int i = 0;
    element_clear(*g);
    for(i = 0; i < NUM_SHARE; i++){
        element_clear(ele_att_hash[i]);
        element_clear(qR[i]);
        element_clear(Cy[i]);
        element_clear(Cy_prime[i]);
    }
    element_clear(qR[i]);
}

int main(void){
    //All names of the variable is either enherited from the paper or from pbc
    //library
    pairing_t pairing;
    element_t g;
    element_t *ele_att_hash;
    element_t *qR;
    element_t *Cy;
    element_t *Cy_prime;
    char s[16384];
    FILE *fp;
    int i = 0;  
    //initialize the pairing
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
    
    //memory allocation starts here
    qR = (element_t *)malloc(sizeof(element_t)*(NUM_SHARE + 1)); 
    ele_att_hash = (element_t *)malloc(sizeof(element_t)*NUM_SHARE);
    Cy = (element_t *)malloc(sizeof(element_t)*NUM_SHARE);
    Cy_prime = (element_t *)malloc(sizeof(element_t)*NUM_SHARE);
    initialize(qR, ele_att_hash, &pairing, &g, Cy, Cy_prime);

    //randomly pick the nodes' values and save them into a file
    fp = fopen("nodes_value_of_tree.txt", "w+");
    if(!fp)
        pbc_die("error opening the tree value file");
    for(i = 0; i < NUM_SHARE + 1; i++){
        element_random(qR[i]);
        fprintf(fp, "qR(%d):", i);
        element_out_str(fp, 10, qR[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    
    element_set(g, ((curve_data_ptr)((a_pairing_data_ptr)
    pairing->data)->Eq->data)->gen);
    /*here is the encryption part*/
    for(i = 0; i < NUM_SHARE; i++){
        element_pow_zn(Cy[i], g, qR[i+1]);    
        element_pow_zn(Cy_prime[i], ele_att_hash[i], qR[i+1]);
    }
    //read out the public key of owner
    //read out the message and multiply it by e(g,g)^(alpha*s)  
    finalize(qR, ele_att_hash, &g, Cy, Cy_prime); 
    return 0;
}
