// Compile with
// gcc -o volume-gmp volume-gmp.c -lgmp
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "time_helper.h"
#include <gmp.h>

#define N 60
#define DIGITS 50


int main() {
    
    unsigned long precision = ceil(DIGITS * log2(10) + 16);
    mpf_set_default_prec(precision);

//	double V[N + 1] = {0.0, 1.0};
    int L[N + 1] = {0, 0, 10, 10, 10, 10,   // 0...5
        10, 10, 10, 10, 10,                 // 6...10
        10, 10, 10, 10, 10,                 // 11...15
        11, 11, 11, 11, 12,                 // 16...20
        12, 12, 12, 13, 13,                 // 21...25
        13, 13, 13, 14, 14,                 // 26...30
        14, 14, 14, 15, 15,                 // 31...35
        15, 15, 15, 16, 16,                 // 36...40
        16, 16, 16, 16, 17,                 // 41...45
        17, 17, 17, 17, 17,                 // 46...50
        18, 18, 18, 18, 18,                 // 51...55
        18, 19, 19, 19, 19};                // 56...60
    
    double V_exact[101] = {0.0,
        2.0000000000,
        3.1415926536,
        4.1887902048,
        4.9348022005,
        5.2637890139,
        5.1677127800,
        4.7247659703,
        4.0587121264,
        3.2985089027,
        2.5501640399,
        1.8841038794,
        1.3352627689,
        9.1062875478e-1,
        5.9926452932e-1,
        3.8144328082e-1,
        2.3533063036e-1,
        1.4098110692e-1,
        8.2145886611e-2,
        4.6621601030e-2,
        2.5806891390e-2,
        1.3949150409e-2,
        7.3704309457e-3,
        3.8106563869e-3,
        1.9295743094e-3,
        9.5772240882e-4,
        4.6630280577e-4,
        2.2287212472e-4,
        1.0463810492e-4,
        4.8287822739e-5,
        2.1915353448e-5,
        9.7871399467e-6,
        4.3030695870e-6,
        1.8634670883e-6,
        7.9520540015e-7,
        3.3452882941e-7,
        1.3878952462e-7,
        5.6808287183e-8,
        2.2948428997e-8,
        9.1522306502e-9,
        3.6047307975e-9,
        1.4025649061e-9,
        5.3926646626e-10,
        2.0494360954e-10,
        7.7007071306e-11,
        2.8615526139e-11,
        1.0518471717e-11,
        3.8254607105e-12,
        1.3768647280e-12,
        4.9053221489e-13,
        1.7302192458e-13,
        6.0433427555e-14,
        2.0906323353e-14,
        7.1644230957e-15,
        2.4325611800e-15,
        8.1846178054e-16,
        2.7293272616e-16,
        9.0220123403e-17,
        2.9567015429e-17,
        9.6079619284e-18,
        3.0962506153e-18,
        9.8964926591e-19,
        3.1377929634e-19,
        9.8700789315e-20,
        3.0805210383e-20,
        9.5408515266e-21,
        2.9326491706e-21,
        8.9473041985e-22,
        2.7097614971e-22,
        8.1474739535e-23,
        2.4322762320e-23,
        7.2101533289e-24,
        2.1225614284e-24,
        6.2058533505e-25,
        1.8022225379e-25,
        5.1990035454e-26,
        1.4899602855e-26,
        4.2423769725e-27,
        1.2002175095e-27,
        3.3741317293e-28,
        9.4264862767e-29,
        2.6173203587e-29,
        7.2229707405e-30,
        1.9813384123e-30,
        5.4027694799e-31,
        1.4646019295e-31,
        3.9472792807e-32,
        1.0577431407e-32,
        2.8183508159e-33,
        7.4674114164e-34,
        1.9675800485e-34,
        5.1559483180e-35,
        1.3437684839e-35,
        3.4834170663e-36,
        8.9820706321e-37,
        2.3038899926e-37,
        5.8787514816e-38,
        1.4923471908e-38,
        3.7691107076e-39,
        9.4714080227e-40,
        2.3682021019e-40};
    
    mpz_t fact[N + 1];
    mpf_t V[N + 1];
    mpf_t error[N + 1];
    for (int i = 0; i <= N; i++) {
        mpf_init(V[i]);
        mpf_init(error[i]);
        mpz_init(fact[i]);
        mpz_fac_ui(fact[i], i);
    }
    mpf_set_ui(V[1], 2.0);
    
    mpz_t P[N + 1][N + 1];
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= i; j++) {
            mpz_init(P[i][j]);
            mpz_tdiv_q(P[i][j], fact[i], fact[j]);
        }
    }
    
    double begin = cpuSecond(), totTime;
    
    printf(" 1 dimension : ");
    gmp_printf("%.*Fe\n", 10, V[1]);
    
    mpz_t m;
    mpz_t cube;
    mpz_t coeff;
    mpz_init(m);
    mpz_init(cube);
    mpz_init(coeff);
    mpf_t numerator, denominator;
    mpf_init(numerator);
    mpf_init(denominator);
	for (int n = 2; n <= N; n++) {
        mpz_set_ui(m, 0);
        mpz_set_ui(cube, L[n] - 1);
        int node[n + 2];
        node[n] = node[n + 1] = 0;
        for (int j = n - 1; j >= 1; j--) {
            mpz_mul_ui(cube, cube, L[n] - 1);
            node[j] = 0;
        }
        int L_sq = (L[n] - 1) * (L[n] - 1);
        while (node[n] < L[n]) {
            int mark = 0;
            int RR = 0;
            for (int j = 1; j <= n; j++) {
                RR += node[j] * node[j];
                if (RR > L_sq) {
                    mark = 1;
                    break;
                }
            }
            if (mark) {
                int k = 1;
                while (node[k] == 0)
                    k++;
                while (k < n && node[k] == node[k + 1])
                    k++;
                node[k]++;
                memset(node, 0, k * sizeof(int));
                continue;
            }
            int edgeCount = 0;
            int k = 1;
            while (k < n && node[k] == node[k + 1]) {
                if (node[k] == 0 || node[k] == L[n] - 1)
                    edgeCount++;
                k++;
            }
            int l = k;
            mpz_set(coeff, P[n][k]);
            int sameCount = 1;
            if (node[k] == 0 || node[k] == L[n] - 1)
                edgeCount++;
            if (k < n) {
                while (++k < n) {
                    if (node[k] == node[k + 1]) {
                        sameCount++;
                        if (k == n - 1)
                            mpz_tdiv_q(coeff, coeff, fact[sameCount]);
                    } else {
                        mpz_tdiv_q(coeff, coeff, fact[sameCount]);
                        sameCount = 1;
                    }
                    if (node[k] == 0 || node[k] == L[n] - 1)
                        edgeCount++;
                }
                if (node[k] == 0 || node[k] == L[n] - 1)
                    edgeCount++;
            }
            mpz_mul_2exp(coeff, coeff, n - edgeCount);
            mpz_add(m, m, coeff);
            node[l]++;
            memset(node, 0, l * sizeof(int));
        }
        mpf_set_z(numerator, m);
        mpf_set_z(denominator, cube);
        mpf_div(V[n], numerator, denominator);
        
        mpf_set_d(error[n], V_exact[n]);
        mpf_div(error[n], V[n], error[n]);
        mpf_sub_ui(error[n], error[n], 1);
        mpf_abs(error[n], error[n]);
        mpf_mul_ui(error[n], error[n], 100);
        
        printf("%2d dimensions: ", n);
        gmp_printf("%.*Fe", 10, V[n]);
        gmp_printf("\t\tError: %.*Ff%%\n", 4, error[n]);
        
        if (n == N) {
            totTime = cpuSecond() - begin;
            printf("Time: %fs\n", totTime);
        } else if (n % 10 == 0)
            printf("Time: %fs\n", cpuSecond() - begin);
	}
    
    printf("Writing to volume.txt ...\n");
    freopen("volume.txt", "w", stdout);
    
    printf(" 1 dimension : ");
    gmp_printf("%.*Fe\n", 10, V[1]);
    for (int n = 2; n <= N; n++) {
        printf("%2d dimensions: ", n);
        gmp_printf("%.*Fe", 10, V[n]);
        gmp_printf("\t\tError: %.*Ff%%\n", 4, error[n]);
    }
    printf("Time: %fs\n", totTime);

	return 0;
}
