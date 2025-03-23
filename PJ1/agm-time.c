#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <time.h>

#define DIGITS 10000
#define M 20

int main(int argc, const char * argv[]) {
    
    float Time[M + 1], average = 0, uncertainty = 0;
    for (int Clock = 1; Clock <= M; Clock++) {
        
        clock_t start = clock();
        
        unsigned long precision = ceil(DIGITS * log2(10));
        
        mpf_t a, b, an, bn, term, sqtwo, S, numerator, denominator, result;
        mpf_set_default_prec(ceil(precision + log2(4.0 * precision) / 4.0));
        mpf_init_set_ui(a, 1L);
        mpf_init(sqtwo);
        mpf_init(b);
        mpf_sqrt_ui(sqtwo, 2L);
        mpf_ui_div(b, 1L, sqtwo);                // a = 1, b = 1/sqrt(2)
        
        unsigned long n = log2(precision * 4);
        mpf_init(S);
        mpf_sub(S, b, a);
        mpf_mul(S, S, S);                    // S = (b - a)^2
        
        mpf_init(an);
        mpf_init(bn);
        mpf_init(term);
        for (unsigned long i = 1; i <= n; i++) {
            mpf_add(an, a, b);
            mpf_div_ui(an, an, 2L);              // an = (a + b) / 2
            mpf_mul(bn, a, b);
            mpf_sqrt(bn, bn);                   // bn = sqrt(a * b)
            mpf_sub(term, bn, an);
            mpf_mul(term, term, term);
            mpf_mul_ui(term, term, 1L << i);     // term = 2^n * (bn - an)^2
            mpf_add(S, S, term);
            mpf_set(a, an);
            mpf_set(b, bn);
        }
        
        mpf_init(numerator);
        mpf_init(denominator);
        mpf_mul(numerator, a, a);
        mpf_mul_ui(numerator, numerator, 4);    // numerator = 4 * an * an
        mpf_ui_sub(denominator, 1L, S);          // denominator = 1 - Sn
        mpf_set_default_prec(precision);
        mpf_init(result);
        mpf_div(result, numerator, denominator);
        
        clock_t stop = clock();
        Time[Clock] = (stop - start) / (float)CLOCKS_PER_SEC;
        
        mpf_clears(a, b, an, bn, term, sqtwo, S, numerator, denominator, result, NULL);
    }
    
    for (int i = 1; i <= M; i++) {
        average += Time[i];
        uncertainty += Time[i] * Time[i];
    }
    average /= M;
    uncertainty /= M;
    uncertainty -= (average * average);
    uncertainty /= (M - 1);
    uncertainty = sqrt(uncertainty);
    printf("Time:(s)\n%f %f\n", average, uncertainty);
    
    return 0;
}
