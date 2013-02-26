#include "stdlib.h"
#include "stdio.h"
#include <pbc/pbc.h>
//#include <pbc/pbc_test.h>
#include <math.h>
#include <gmp.h>
#include "grs.h"
#include "string.h"
#include "stdbool.h"
#include "message_handle.h"
#include <string.h>

using namespace std;

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

void initialize(element_t *qR, 
                element_t *ele_att_hash, 
                pairing_t *pairing,
                element_t *g,
                element_t *h,
                element_t *C,
                element_t *C_tilde,
                element_t *Cy,
                element_t *Cy_prime,
                mpz_t *message_mpz,
                element_t *message_ele){
   //initialize the attribute string 
    int i = 0;
    char *attr[NUM_SHARE];
    strcpy(attr[0], "http://www.servername.com/pathtofile");
    strcpy(attr[1], "Julia Zhu");
    strcpy(attr[2], "4");
    strcpy(attr[3], "r");
    strcpy(attr[4], "2013/01/10/13/20");
    strcpy(attr[5], "julia@facebook");
    strcpy(attr[6], "printer");
    strcpy(attr[7], "check string 1");
    strcpy(attr[8], "check string 2");
    element_init_G1(*g, *pairing);
    element_init_G1(*h, *pairing);
    element_init_G1(*C, *pairing);
    element_init_GT(*C_tilde, *pairing);
    mpz_init(*message_mpz);
    element_init_GT(*message_ele, *pairing);
    //initialize the attributes' hash value
    for(i = 0; i < NUM_SHARE; i++){
        element_init_G1(ele_att_hash[i], *pairing);
        element_from_hash(ele_att_hash[i], (void *)attr[i], 256);
        element_init_Zr(qR[i], *pairing);
        element_init_G1(Cy[i], *pairing);
        element_init_G1(Cy_prime[i], *pairing);
    }
    element_init_Zr(qR[i], *pairing);
}

void finalize(element_t *qR,
                element_t *ele_att_hash,
                element_t *g,
                element_t *h,
                element_t *C,
                element_t *C_tilde,
                element_t *Cy,
                element_t *Cy_prime,
                mpz_t *message_mpz,
                element_t *message_ele){
    int i = 0;
    element_clear(*g);
    element_clear(*h);
    element_clear(*C);
    element_clear(*C_tilde);
    mpz_clear(*message_mpz);
    element_clear(*message_ele);
    for(i = 0; i < NUM_SHARE; i++){
        element_clear(ele_att_hash[i]);
        element_clear(qR[i]);
        element_clear(Cy[i]);
        element_clear(Cy_prime[i]);
    }
    element_clear(qR[i]);
}

int encryption(void){
    //All names of the variable is either enherited from the paper or from pbc
    //library
    pairing_t pairing;
    element_t g;
    element_t h;
    element_t *ele_att_hash;
    element_t *qR;
    element_t *C;
    element_t *C_tilde;
    element_t *Cy;
    element_t *Cy_prime;
    element_t *message_ele;
    char s[16384];
    FILE *fp, *fp_write;
    int i = 0; 
    int num_read = 0;
    char *str_read, *str1;
    char *m;
    mpz_t message_mpz;
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
    if(!pairing_is_symmetric(pairing)) printf("error\n");//pbc_die("pairing is not symmetric\n");
    
    //memory allocation starts here
    qR = (element_t *)malloc(sizeof(element_t)*(NUM_SHARE + 1)); 
    ele_att_hash = (element_t *)malloc(sizeof(element_t)*NUM_SHARE);
    C = (element_t *)malloc(sizeof(element_t));
    C_tilde = (element_t *)malloc(sizeof(element_t));
    Cy = (element_t *)malloc(sizeof(element_t)*NUM_SHARE);
    Cy_prime = (element_t *)malloc(sizeof(element_t)*NUM_SHARE);
    str_read = (char *)malloc(sizeof(char)*500);
    str1 = (char *)malloc(sizeof(char)*3);
    m = (char *)malloc(sizeof(char)*512);
    message_ele = (element_t *)malloc(sizeof(element_t));
    //initialize the elements to the coresponding type
    initialize(qR, ele_att_hash, &pairing, &g, &h, C, C_tilde, Cy, Cy_prime,
    &message_mpz, message_ele);

    //randomly pick the nodes' values qR[i] and save them into a file
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
    
 //   curve_data_ptr p3 = (curve_data_ptr) (*(*(a_pairing_data_ptr) (*pairing).data).Eq).data;

    element_set(g, (*(curve_data_ptr) (*(*(a_pairing_data_ptr) (*pairing).data).Eq).data).gen);
    //here is the encryption part
    //compute Cs and C_primes
    for(i = 0; i < NUM_SHARE; i++){
        element_pow_zn(Cy[i], g, qR[i+1]);   
        element_pow_zn(Cy_prime[i], ele_att_hash[i], qR[i+1]);
    }
    //read out the public key h of authorizer
    fp = fopen("../public/authorizer_public_keys.txt", "r");
    if(!fp)
        pbc_die("open authorizer public key file failed");
    else{
        while(true){
            fgets(str_read, 500, fp);
            strncpy(str1, str_read, 2);
            if(strstr(str_read, "h:"))
                break;
        }
        str_read += 2;
        element_set_str(h, str_read, 10);
        fclose(fp);
    }
    //compute C
    element_pow_zn(*C, h, qR[0]);
    //read out the owner public key e(g, g)^alpha 
    fp = fopen("../public/owner_public_key.txt", "r");
    if(!fp)
        pbc_die("open owner public key file failed");
    else{
        while(true){
            fgets(str_read, 500, fp);
            strncpy(str1, str_read, 2);
            if(strstr(str_read, "public key:"))
                break;
        }
        str_read += 11;
        element_set_str(*C_tilde, str_read, 10);
        //printf("%s\n", str_read);
        fclose(fp);
    }
    element_pow_zn(*C_tilde, *C_tilde, qR[0]);//e(g, g)^(alpha*s)
    
    //read out the message and  
    fp = fopen("password.txt", "r");
    if(!fp)
        pbc_die("fail to open the message file");
    else{
        fp_write = fopen("encryption_header.txt", "w+");
        if(!fp_write)
            pbc_die("creat encryption file failed");
        fprintf(fp_write, "C_tilde:");
        while(true){
            num_read = fread(str_read, 1, 100, fp);
            str_read[num_read] = '\0';
           // printf("%s%d\n", str_read, num_read);
           //we need to initialize m_mpz to 0 in order to start a another round
           //of encryption.
           mpz_init(message_mpz);
           messageToValue(str_read, message_mpz, m);
           
           strcpy(str_read, "[");
           strcat(str_read, m);
           strcat(str_read, ",0]");
           element_set_str(*message_ele, str_read, 10);
           //multiply it by e(g,g)^(alpha*s)
           element_mul(*C_tilde, *C_tilde, *message_ele);
           //write the encrypted data to file
           element_out_str(fp_write, 10, *C_tilde);
           if(feof(fp)){
                break;
            }
     /*       if(ferror(fp))
                pbc_die("reading error occurs");
        */
        }
        //write the Cy and Cy_prime into the file
        // i starts from 1, because qR[0] is secret s that we want to protect, qR[1], qR[2], ...,
        // qR[NUM_SHARE] are the secret shares
        for(i = 1; i < NUM_SHARE; i++){
            fprintf(fp_write,"\nCy[%d]:", i);
            element_out_str(fp_write, 10, Cy[i]);
            fprintf(fp_write, "\nCy_prime[%d]:", i);
            element_out_str(fp_write, 10, Cy_prime[i]);
        }
        fclose(fp);
        fclose(fp_write);
    }
    finalize(qR, ele_att_hash, &g, &h, C, C_tilde, Cy, Cy_prime, &message_mpz,
    message_ele); 
    return 0;
}
