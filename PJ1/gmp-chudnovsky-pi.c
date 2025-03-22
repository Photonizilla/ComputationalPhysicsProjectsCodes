// The GNU Multiple Precision Arithmetic Library (gmp) required
// https://gmplib.org/
//
// # For Debian/Ubuntu
// sudo apt-get install libgmp3-dev
//
// # For Fedora
// sudo dnf install gmp-devel
//
// # For Arch Linux
// sudo pacman -S gmp
//
// Compile with
// gcc -o gmp-chudnovksy-pi gmp-chudnovsky-pi.c -lgmp
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>

#define DIGITS 100000000
#define DIGITS_PER_ITER  14.1816474627254776555
#define A 13591409
#define B 545140134
#define C 640320
#define D 53360

mpz_t* PStack;
mpz_t* QStack;
mpz_t* RStack;
unsigned long top;

void binarySplitting(unsigned long a, unsigned long b, int level) {
    if (b == a + 1) {
        mpz_set_ui(PStack[top], (2 * a - 1));
        mpz_mul_ui(PStack[top], PStack[top], (6 * a - 5));
        mpz_mul_ui(PStack[top], PStack[top], (6 * a - 1));
        mpz_neg(PStack[top], PStack[top]);                          // P(a, a+1) = -(6a - 1)(2a - 1)(6a - 5)
        mpz_set_ui(QStack[top], a);
        mpz_mul_ui(QStack[top], QStack[top], a);
        mpz_mul_ui(QStack[top], QStack[top], a);
        mpz_mul_ui(QStack[top], QStack[top], (C / 24) * (C / 24));
        mpz_mul_ui(QStack[top], QStack[top], C * 24);               // Q(a, a+1) = C^3 * a^3 / 24
        mpz_set_ui(RStack[top], a);
        mpz_mul_ui(RStack[top], RStack[top], B);
        mpz_add_ui(RStack[top], RStack[top], A);
        mpz_mul(RStack[top], RStack[top], PStack[top]);             // R(a, a+1) = P(a, a+1) * (A + Ba)
        return;
    }
    unsigned long mid = a + (b-a) * 0.5224;
    binarySplitting(a, mid, level + 1);
    top++;
    binarySplitting(mid, b, level + 1);
    top--;
    
    mpz_mul(RStack[top], RStack[top], QStack[top + 1]);
    mpz_mul(RStack[top + 1], RStack[top + 1], PStack[top]);
    mpz_add(RStack[top], RStack[top], RStack[top + 1]);             // R(a, b) = R(a, mid) Q(mid, b) + P(a, mid) R(mid, b)
    mpz_mul(PStack[top], PStack[top], PStack[top + 1]);             // P(a, b) = P(a, mid) P(mid, b)
    mpz_mul(QStack[top], QStack[top], QStack[top + 1]);             // Q(a, b) = Q(a, mid) Q(mid, b)
}

int main(int argc, const char * argv[]) {
    
    clock_t start = clock();
    
    unsigned long precision = ceil(DIGITS * log2(10) + 16);
    
    mpf_t result, res_numerator, res_denominator, Q, R;
    
    unsigned long terms = ceil((DIGITS + 1) / DIGITS_PER_ITER);
    unsigned long depth = 1;
    while ((1L << depth) < terms)
        depth++;
    depth++;
    
    PStack = (mpz_t*)malloc(sizeof(mpz_t) * depth);
    QStack = (mpz_t*)malloc(sizeof(mpz_t) * depth);
    RStack = (mpz_t*)malloc(sizeof(mpz_t) * depth);
    
    for (unsigned long i = 0; i < depth; i++) {
        mpz_init(PStack[i]);
        mpz_init(QStack[i]);
        mpz_init(RStack[i]);
    }
    
    if (terms <= 1) {
        mpz_set_ui(PStack[top], 1);
        mpz_set_ui(QStack[top], 1);
        mpz_set_ui(RStack[top], 0);
    } else {
        binarySplitting(1, terms, 0);
    }
    
    mpf_set_default_prec(precision);
    mpf_init(Q);
    mpf_set_z(Q, QStack[top]);
    mpf_init(R);
    mpf_set_z(R, RStack[top]);
    mpf_init(res_numerator);
    mpf_sqrt_ui(res_numerator, C);
    mpf_mul_ui(res_numerator, res_numerator, D);
    mpf_mul(res_numerator, res_numerator, Q);
    mpf_init_set_ui(res_denominator, A);
    mpf_mul(res_denominator, res_denominator, Q);
    mpf_add(res_denominator, res_denominator, R);
    mpf_init(result);
    mpf_div(result, res_numerator, res_denominator);
    
    clock_t stop = clock();
    printf("Time:(s)\n%f\n", (stop - start) / (float)CLOCKS_PER_SEC);
    
    freopen("gmp-chudnovsky-pi.txt", "w", stdout);
    gmp_printf("%.*Ff\n", DIGITS, result);
    
    for (unsigned long i = 0; i < depth; i++) {
        mpz_clear(PStack[i]);
        mpz_clear(QStack[i]);
        mpz_clear(RStack[i]);
    }
    mpf_clears(result, res_numerator, res_denominator, Q, R, NULL);
    
    free(PStack);
    free(QStack);
    free(RStack);
    
    return 0;
}


//#include <stdio.h>
//#include <math.h>
//#include <gmp.h>
//#include <time.h>
//
//#define DIGITS 10000
//#define DIGITS_PER_ITER  14.1816474627254776555
//#define A 13591409
//#define B 545140134
//#define C 640320
//#define D 53360
//
//int main(int argc, const char * argv[]) {
//    
//    clock_t start = clock();
//    
//    int precision = ceil(DIGITS * log2(10) + 100);      // ?
//    
//    mpf_t result, coefficient, term, numerator, denominator;
//    mpz_t numerator_add, numerator_fac, denominator_mul, denominator_fac, denominator_cub;
//    mpf_set_default_prec(precision);
//    mpf_init_set_ui(result, A);
//    mpf_init(numerator);
//    mpz_init_set_ui(numerator_add, A);
//    mpz_init(numerator_fac);
//    mpf_init(denominator);
//    mpz_init_set_ui(denominator_mul, 1);
//    mpz_init(denominator_fac);
//    mpz_init(denominator_cub);
//    mpf_init(term);
//    int n = ceil((DIGITS + 1) / DIGITS_PER_ITER);      // ?
//    for (int i = 1; i <= n; i++) {
//        int j = 3 * i;
//        mpz_add_ui(numerator_add, numerator_add, B);                // (B * k + A)
//        mpz_fac_ui(numerator_fac, j << 1);                          // (6 * k)!
//        mpz_mul(numerator_fac, numerator_fac, numerator_add);
//        mpf_set_z(numerator, numerator_fac);                        // numerator = (6 * k)! * (B * k + A)
//        mpz_mul_ui(denominator_mul, denominator_mul, C);            // C^k
//        mpz_mul_ui(denominator_mul, denominator_mul, i);            // k!
//        mpz_pow_ui(denominator_cub, denominator_mul, 3);            // (k!)^3 * C^(3 * k)
//        mpz_fac_ui(denominator_fac, j);                             // (3 * k)!
//        mpz_mul(denominator_fac, denominator_fac, denominator_cub);
//        mpf_set_z(denominator, denominator_fac);                    // denominator = (3 * k)! * (k!)^3 * C^(3 * k)
//        mpf_div(term, numerator, denominator);
//        if (i % 2 == 0)                                                  // (-1)^k
//            mpf_add(result, result, term);
//        else
//            mpf_sub(result, result, term);
//    }
//    
//    mpf_init(coefficient);
//    mpf_sqrt_ui(coefficient, C);
//    mpf_mul_ui(coefficient, coefficient, D);
//    mpf_div(result, coefficient, result);
//    
//    clock_t stop = clock();
//    printf("Time: %fs\n", (stop - start) / (float)CLOCKS_PER_SEC);
//    
//    freopen("gmp-chudnovsky-pi.txt", "w", stdout);
//    gmp_printf("%.*Ff\n", DIGITS, result);
//    
//    mpf_clears(result, coefficient, term, numerator, denominator, NULL);
//    mpz_clears(numerator_add, numerator_fac, denominator_mul, denominator_fac, denominator_cub, NULL);
//    
//    return 0;
//}
