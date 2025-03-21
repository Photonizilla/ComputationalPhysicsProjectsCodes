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
// gcc -o gmp-agm-pi gmp-agm-pi.c -lgmp
//

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <time.h>

#define DIGITS 100000

int main(int argc, const char * argv[]) {
    
    clock_t start = clock();
    
    int precision = ceil(DIGITS * log2(10) + 100);
    
    mpf_t a, b, an, bn, term, sqtwo, S, numerator, denominator, result;
    mpf_set_default_prec(ceil(precision + log2(4.0 * precision) / 4.0));
    mpf_init_set_ui(a, 1);
    mpf_init(sqtwo);
    mpf_init(b);
    mpf_sqrt_ui(sqtwo, 2);
    mpf_ui_div(b, 1, sqtwo);                // a = 1, b = 1/sqrt(2)
    
    int n = log2(precision * 4);
    mpf_init(S);
    mpf_sub(S, b, a);
    mpf_pow_ui(S, S, 2);                    // S = (b - a)^2
    
    mpf_init(an);
    mpf_init(bn);
    mpf_init(term);
    for (int i = 1; i <= n; i++) {
        mpf_add(an, a, b);
        mpf_div_ui(an, an, 2);              // an = (a + b) / 2
        mpf_mul(bn, a, b);
        mpf_sqrt(bn, bn);                   // bn = sqrt(a * b)
        mpf_sub(term, bn, an);
        mpf_pow_ui(term, term, 2);
        mpf_mul_ui(term, term, 1 << i);     // term = 2^n * (bn - an)^2
        mpf_add(S, S, term);
        mpf_set(a, an);
        mpf_set(b, bn);
    }
    
    mpf_init(numerator);
    mpf_init(denominator);
    mpf_mul(numerator, a, a);
    mpf_mul_ui(numerator, numerator, 4);    // numerator = 4 * an * an
    mpf_ui_sub(denominator, 1, S);          // denominator = 1 - Sn
    mpf_set_default_prec(precision);
    mpf_init(result);
    mpf_div(result, numerator, denominator);
    
    clock_t stop = clock();
    printf("Time: %fs\n", (stop - start) / (float)CLOCKS_PER_SEC);
    
    freopen("gmp-agm-pi.txt", "w", stdout);
    gmp_printf("%.*Ff\n", DIGITS, result);
    
    mpf_clears(a, b, an, bn, term, sqtwo, S, numerator, denominator, result, NULL);
    
    return 0;
}
