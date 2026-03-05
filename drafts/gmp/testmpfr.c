#include <stdio.h> 
#include <mpfr.h>

int main() {
    printf("mpfr_int = %lu\n", sizeof(mpfr_int));
    printf("mpfr_uint = %lu\n", sizeof(mpfr_uint));
    printf("mpfr_long = %lu\n", sizeof(mpfr_long));
    printf("mpfr_ulong = %lu\n", sizeof(mpfr_ulong));
    printf("mpfr_size_t = %lu\n", sizeof(mpfr_size_t));
    printf("mpfr_prec_t = %lu\n", sizeof(mpfr_prec_t));
    printf("mpfr_sign_t = %lu\n", sizeof(mpfr_sign_t));
    printf("mpfr_exp_t = %lu\n", sizeof(mpfr_exp_t));
    printf("mpfr_uexp_t = %lu\n", sizeof(mpfr_uexp_t));
    printf("mp_limb_t = %lu\n", sizeof(mp_limb_t));
    printf("mpfr_t = %lu\n", sizeof(mpfr_t));
    printf("mpfr_t* = %lu\n", sizeof(mpfr_t*));

/*
    typedef int             mpfr_int;
typedef unsigned int    mpfr_uint;
typedef long            mpfr_long;
typedef unsigned long   mpfr_ulong;
typedef size_t          mpfr_size_t;
typedef long  mpfr_prec_t;
typedef unsigned long  mpfr_uprec_t;
typedef int          mpfr_sign_t;
typedef short mpfr_exp_t;
typedef unsigned short mpfr_uexp_t;
  mp_limb_t   *_mpfr_d;
typedef __mpfr_struct mpfr_t[1];
*/
}
