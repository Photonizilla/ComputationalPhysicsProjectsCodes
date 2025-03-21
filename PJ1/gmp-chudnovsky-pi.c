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
#include <math.h>
#include <gmp.h>
#include <time.h>

#define DIGITS 10000
#define DIGITS_PER_ITER  14.1816474627254776555
#define A 13591409
#define B 545140134
#define C 640320
#define D 53360

int main(int argc, const char * argv[]) {
    
    clock_t start = clock();
    
    int precision = ceil(DIGITS * log2(10) + 100);      // ?
    
    mpf_t result, coefficient, term, numerator, denominator;
    mpz_t numerator_add, numerator_fac, denominator_mul, denominator_fac, denominator_cub;
    mpf_set_default_prec(precision);
    mpf_init_set_ui(result, A);
    mpf_init(numerator);
    mpz_init_set_ui(numerator_add, A);
    mpz_init(numerator_fac);
    mpf_init(denominator);
    mpz_init_set_ui(denominator_mul, 1);
    mpz_init(denominator_fac);
    mpz_init(denominator_cub);
    mpf_init(term);
    int n = ceil((DIGITS + 1) / DIGITS_PER_ITER);      // ?
    for (int i = 1; i <= n; i++) {
        int j = 3 * i;
        mpz_add_ui(numerator_add, numerator_add, B);                // (B * k + A)
        mpz_fac_ui(numerator_fac, j << 1);                          // (6 * k)!
        mpz_mul(numerator_fac, numerator_fac, numerator_add);
        mpf_set_z(numerator, numerator_fac);                        // numerator = (6 * k)! * (B * k + A)
        mpz_mul_ui(denominator_mul, denominator_mul, C);            // C^k
        mpz_mul_ui(denominator_mul, denominator_mul, i);            // k!
        mpz_pow_ui(denominator_cub, denominator_mul, 3);            // (k!)^3 * C^(3 * k)
        mpz_fac_ui(denominator_fac, j);                             // (3 * k)!
        mpz_mul(denominator_fac, denominator_fac, denominator_cub);
        mpf_set_z(denominator, denominator_fac);                    // denominator = (3 * k)! * (k!)^3 * C^(3 * k)
        mpf_div(term, numerator, denominator);
        if (i % 2 == 0)                                                  // (-1)^k
            mpf_add(result, result, term);
        else
            mpf_sub(result, result, term);
    }
    
    mpf_init(coefficient);
    mpf_sqrt_ui(coefficient, C);
    mpf_mul_ui(coefficient, coefficient, D);
    mpf_div(result, coefficient, result);
    
    clock_t stop = clock();
    printf("Time: %fs\n", (stop - start) / (float)CLOCKS_PER_SEC);
    
    freopen("gmp-chudnovsky-pi.txt", "w", stdout);
    gmp_printf("%.*Ff\n", DIGITS, result);
    
    mpf_clears(result, coefficient, term, numerator, denominator, NULL);
    mpz_clears(numerator_add, numerator_fac, denominator_mul, denominator_fac, denominator_cub, NULL);
    
    return 0;
}
