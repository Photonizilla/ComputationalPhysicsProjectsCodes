// Compile with
/*
g++ -c sobol_generator.cpp -o sobol_generator.o
g++ -c volume_sobol.cpp -o volume_sobol.o
g++ sobol_generator.o volume_sobol.o -o volume_sobol -lgmp -std=c++11
*/
//

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "time_helper.h"
#include <gmp.h>
#include "sobol_generator.hpp"

#define DIGITS 10
#define N 200
#define GENERATE_SIZE_IN_G 4        // <= 6
#define GIGABYTE_SIZE 134217728

const double V[N + 1] = {0.0,
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
        2.3682021019e-40,
        5.8921397749e-41,
        1.4588090834e-41,
        3.5943112681e-42,
        8.8134305759e-43,
        2.1508308333e-43,
        5.2241903302e-44,
        1.2629970738e-44,
        3.0393107338e-45,
        7.2804079424e-46,
        1.7360502315e-46,
        4.1210947940e-47,
        9.7392190242e-48,
        2.2914692265e-48,
        5.3678349014e-49,
        1.2519761544e-49,
        2.9075087400e-50,
        6.7234172464e-51,
        1.5481708640e-51,
        3.5499560048e-52,
        8.1062036883e-53,
        1.8433910257e-53,
        4.1748180255e-54,
        9.4165588681e-55,
        2.1154157482e-55,
        4.7332787459e-56,
        1.0548848530e-56,
        2.3417375962e-57,
        5.1781539134e-58,
        1.1405869192e-58,
        2.5027154298e-59,
        5.4706251696e-60,
        1.1912897588e-60,
        2.5844324577e-61,
        5.5858912755e-62,
        1.2028494848e-62,
        2.5806757345e-63,
        5.5165884744e-64,
        1.1749901346e-64,
        2.4936509099e-65,
        5.2733433928e-66,
        1.1112106921e-66,
        2.3333375863e-67,
        4.8824774082e-68,
        1.0181105860e-68,
        2.1156903665e-69,
        4.3814914214e-70,
        9.0430439627e-71,
        1.8601163867e-71,
        3.8133638227e-72,
        7.7916373002e-73,
        1.5867597048e-73,
        3.2208092766e-74,
        6.5162779496e-75,
        1.3140871119e-75,
        2.6414827013e-76,
        5.2927261757e-77,
        1.0571289999e-77,
        2.1047581862e-78,
        4.1774449056e-79,
        8.2653660692e-80,
        1.6302894691e-80,
        3.2057300398e-81,
        6.2843011280e-82,
        1.2281826759e-82,
        2.3930562735e-83,
        4.6487345444e-84,
        9.0036024054e-85,
        1.7386226539e-85,
        3.3474143400e-86,
        6.4259343024e-87,
        1.2299663508e-87,
        2.3474032554e-88,
        4.4671135858e-89,
        8.4765342785e-90,
        1.6038687113e-90,
        3.0261156611e-91,
        5.6934487692e-92,
        1.0681823292e-92,
        1.9984912656e-93,
        3.7286597312e-94,
        6.9375088157e-95,
        1.2872450570e-95,
        2.3819482765e-96,
        4.3956517550e-97,
        8.0898499530e-98,
        1.4848760496e-98,
        2.7181832279e-99,
        4.9626337116e-100,
        9.0364279999e-101,
        1.6411130117e-101,
        2.9726466826e-102,
        5.3705297720e-103,
        9.6775595647e-104,
        1.7393831833e-104,
        3.1182512854e-105,
        5.5759524799e-106,
        9.9454571882e-107,
        1.7694314493e-107,
        3.1401583155e-108,
        5.5588328420e-109};

int main(int argc, const char* argv[]) {
    
//    Xorshift32 start(time(NULL));
    
	double begin = cpuSecond();

    int dim, logCount;
    double generate_size_in_G;
    if (argc < 4) {
        printf("Usage: command <dimensions> <log_2(# of samples)> <memory size in GB>\n");
        return -1;
    }
    dim = atoi(argv[1]);
    unsigned long precision = ceil((DIGITS - log10(V[dim])) * log2(10) + 16);
    mpf_set_default_prec(precision);
    logCount = atoi(argv[2]);
    generate_size_in_G = atof(argv[3]);
    if (dim < 2) {
        printf("Dimensions should not be smaller than 2.\n");
        return 1;
    }
    if (logCount < 0) {
        printf("log_2(# of samples) should be a natural number.\n");
        return 2;
    }
    if (generate_size_in_G <= 0.0) {
        printf("memory size must be greater than 0.\n");
        return 3;
    }
    int block_size = generate_size_in_G * GIGABYTE_SIZE / dim;
    long long shift = block_size;
    int blockIdx = 0;
    double* samples = new double[(long long)(generate_size_in_G * GIGABYTE_SIZE)];
    
    if (logCount <= 30) {
        int targetCount = 1 << logCount;
        int inBall = 0;
        int total = 0;
        double ratio;
        if (targetCount <= block_size) {
            total = targetCount;
            printf("Generating %d samples of dimension %d... ",
                   total, dim);
            i8_sobol_generate(samples, dim, total, total);
            printf("Analysing...\n");
            for (int i = 0; i < total; i++) {
                double RR = 0.0;
                int mark = 1;
                for (int d = 0; d < dim; d++) {
                    RR += samples[i * dim + d] * samples[i * dim + d];
                    if (RR > 1.0) {
                        mark = 0;
                        break;
                    }
                }
                if (mark) {
                    inBall++;
                }
            }
            ratio = (double)inBall / total;
        } else {
            int block_inBall = 0;
            while (total < targetCount) {
                printf("Block %d: Generating %d samples of dimension %d... ",
                       blockIdx++, block_size, dim);
//                shift = (shift + start.next() % block_size) % 0x3fffffff;
                shift = (shift + block_size) % LLONG_MAX;
                i8_sobol_generate(samples, dim, block_size,
                                  block_size + shift);
                printf("Analysing...\n");
                block_inBall = 0;
                for (int i = 0; i < block_size; i++) {
                    double RR = 0.0;
                    int mark = 1;
                    for (int d = 0; d < dim; d++) {
                        RR += samples[i * dim + d] * samples[i * dim + d];
                        if (RR > 1.0) {
                            mark = 0;
                            break;
                        }
                    }
                    if (mark) {
                        block_inBall++;
                    }
                }
                inBall += block_inBall;
                total += block_size;
            }
            ratio = (double)inBall / total;
        }
        delete[] samples;
        printf("%d samples in total.\n", total);
        if (dim <= 30) {
            double result = ratio * (1 << dim);
            printf("V = %g\n", result);
            printf("Error: %.4f%%\n", 100.0 * fabs(result - V[dim]) / V[dim]);
        } else {
            printf("Using gmp library.\n");
            mpz_t cubeVolume;
            mpz_init_set_ui(cubeVolume, 1);
            mpz_mul_2exp(cubeVolume, cubeVolume, dim);
            mpf_t result;
            mpf_init(result);
            mpf_set_z(result, cubeVolume);
            mpf_t ratio_mpf;
            mpf_init_set_d(ratio_mpf, ratio);
            mpf_mul(result, result, ratio_mpf);
            mpf_t error;
            mpf_init(error);
            mpf_set_d(error, V[dim]);
            mpf_div(error, result, error);
            mpf_sub_ui(error, error, 1);
            mpf_abs(error, error);
            mpf_mul_ui(error, error, 100);
            gmp_printf("V = %.*Fe\n", 10, result);
            gmp_printf("Error: %.*Ff%%\n", 4, error);
            mpz_clear(cubeVolume);
            mpf_clears(result, ratio_mpf, error, NULL);
        }
        printf("Time: %fs\n", cpuSecond() - begin);
        return 0;
    }
    
    printf("Using gmp library.\n");
    mpz_t targetCount;
    mpz_t total;
    mpz_t inBall;
    mpz_init_set_ui(targetCount, 1);
    mpz_init_set_ui(total, 0);
    mpz_init_set_ui(inBall, 0);
    mpz_mul_2exp(targetCount, targetCount, logCount);
    mpf_t numerator;
    mpf_t denominator;
    mpf_t result;
    mpf_init(numerator);
    mpf_init(denominator);
    mpf_init(result);
    
    int block_inBall = 0;
    while (mpz_cmp(total, targetCount) < 0) {
        printf("Block %d: Generating %d samples of dimension %d... ",
               blockIdx++, block_size, dim);
        
//        shift = (shift + start.next() % block_size) % 0x3fffffff;
        shift = (shift + block_size) % LLONG_MAX;
        i8_sobol_generate(samples, dim, block_size,
                          block_size + shift);
        printf("Analysing...\n");
        block_inBall = 0;
        for (int i = 0; i < block_size; i++) {
            double RR = 0.0;
            int mark = 1;
            for (int d = 0; d < dim; d++) {
                RR += samples[i * dim + d] * samples[i * dim + d];
                if (RR > 1.0) {
                    mark = 0;
                    break;
                }
            }
            if (mark) {
                block_inBall++;
            }
        }
        mpz_add_ui(inBall, inBall, block_inBall);
        mpz_add_ui(total, total, block_size);
    }
    delete[] samples;
    gmp_printf ("%Zd samples in total.\n", total);
    mpf_set_z(numerator, inBall);
    mpf_set_z(denominator, total);
    mpf_div(result, numerator, denominator);
    mpz_t cubeVolume;
    mpz_init_set_ui(cubeVolume, 1);
    mpz_mul_2exp(cubeVolume, cubeVolume, dim);
    mpf_set_z(numerator, cubeVolume);
    mpf_mul(result, result, numerator);
    mpf_t error;
    mpf_init(error);
    mpf_set_d(error, V[dim]);
    mpf_div(error, result, error);
    mpf_sub_ui(error, error, 1);
    mpf_abs(error, error);
    mpf_mul_ui(error, error, 100);
    gmp_printf("V = %.*Fe\n", 10, result);
    gmp_printf("Error: %.*Ff%%\n", 4, error);
	printf("Time: %fs\n", cpuSecond() - begin);
    
    mpz_clears(targetCount, total, inBall, cubeVolume, NULL);
    mpf_clears(result, numerator, denominator, error, NULL);

	return 0;
}

