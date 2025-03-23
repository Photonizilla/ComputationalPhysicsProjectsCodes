#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>

#define DIGITS 10000
#define DIGITS_PER_ITER  14.1816474627254776555
#define A 13591409
#define B 545140134
#define C 640320
#define D 53360

#define M 20

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
    
    float Time[M + 1], average = 0, uncertainty = 0;
    for (int Clock = 1; Clock <= M; Clock++) {
        
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
        Time[Clock] = (stop - start) / (float)CLOCKS_PER_SEC;
        
        for (unsigned long i = 0; i < depth; i++) {
            mpz_clear(PStack[i]);
            mpz_clear(QStack[i]);
            mpz_clear(RStack[i]);
        }
        mpf_clears(result, res_numerator, res_denominator, Q, R, NULL);
        
        free(PStack);
        free(QStack);
        free(RStack);
        top = 0;
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
