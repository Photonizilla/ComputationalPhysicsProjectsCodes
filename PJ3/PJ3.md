# Volume of n-dimensional sphere using integration

## Result

Using the **trapezoidal** method on an $n$-dimensional space, we can integrate the volume of spheres of dimensions 1 to 60 within 35 seconds.

| Dimensions | Time                      | Error             |
| ---------- | ------------------------- | ----------------- |
| 1 ~ 10     | $ 0.000732 \pm 0.000067 $ | $ \leq 3.2762\% $ |
| 1 ~ 20     | $ 0.01252 \pm 0.00070 $   | $ \leq 4.7977\% $ |
| 1 ~ 30     | $ 0.1575 \pm 0.0013 $     | $ \leq 4.9336\% $ |
| 1 ~ 40     | $ 1.2650 \pm 0.0014 $     | $ \leq 4.9564\% $ |
| 1 ~ 50     | $ 6.5638 \pm 0.0017 $     | $ \leq 4.9649\% $ |
| 1 ~ 60     | $ 34.138 \pm 0.049 $      | $ \leq 4.9649\% $ |

These data are obtained by averaging 20 calculations.
Furthermore, we can calculate the volume of a 150-dimensional sphere in about 1 minute within the error of 10%.

| Dimensions | Time           | Volume (numerical)       | Error        |
| ---------- | -------------- | ------------------------ | ------------ |
| 150        | $ 64.975598 $s | $8.5446 \times 10^{-73}$ | $ 9.6639\% $ |

## Method

A point in $n$-dimensional space is given by $ (x_1, x_2, \ldots, x_n) $. We calculate the volume of an $n$-dimensional sphere of radius $ R $ by integrating
$$  f(x_1, x_2, \ldots, x_n) =
    \begin{cases}
        1, & \left( \sum_{i = 1}^n x_i^2 \right) \leq R^2,\\
        0, & \text{otherwise}
    \end{cases}
$$ over the region $ -R \leq x_i \leq R,\ i = 1, 2, \ldots, n $.

We use the trapezoidal method. If we sample $ L $ points on each dimension, the naive algorithm requires $ L^n $ iterations.

Remark
: We considered the Romberg method, despite that it converges very fast, it involves the calculation of a large amount of `float`s or `double`s, which actually slows it down.

### Simplifications

1. Using the fact that the sphere is an even function for each coordinate, thus is symmetrical with respect to each quadrant, octant, etc., that is, using the reflective symmetry:
$$  f(x_i) = f(-x_i), \quad \forall i = 1, 2, \ldots, n.
$$ **Now, we only have to calculate the volume of the region $ 0 \leq x_i \leq R $ and multiply the result by $ 2^n $.**
We again sample $ L $ points evenly on $ [0, R] $, the first point would be $0$, the last would be $R$.
To compute the volume of a sphere with radius $ R = 1 $, we can first set $ R = L - 1 $, and divide the result by $ (L - 1)^n $. With $ R $ set to $ (L - 1) $, we assign volume $ 1 $ to points inside the $n$-dimensional cube, and assign volume $ 1/2 $ to points on surface, then assign volume $ 1/4 $ to points on $(n-2)$-dimensional surface, and volume $ 1/2^m $ to points on $(n-m)$-dimensional surface. In this way, total volume of the n-dimensional cube will be $ (L - 1)^n $.

2. Consider a $1 / 8$ ball in a cube of $ R \times R \times R $. We see that many verticle lines inside the cube can be drawn without touching the ball. This implies that when we calculate a 3-dimensional sphere, many 1-dimensional spheres in the integration are unnecessary in the iterations. By the same token, when we calculate an $n$-dimensional sphere, many spheres of dimension $n - 2, n - 3, \ldots$ are unnecessary in the iterations. This tells us that we can skip many regions in our iterations.
Mathematically, if point $(a, b, c, 0, 0)$ is outside of the sphere (i.e., $ a^2 + b^2 + c^2 > R^5 $), we know that
$$  \begin{aligned}
        &(a, b, c, 0, 1), (a, b, c, 0, 2), \ldots\\
        &(a, b, c, 1, \#), (a, b, c, 2, \#), \ldots\\
        &(a, b, c + 1, \#, \#), (a, b, c + 2, \#, \#), \ldots
    \end{aligned}
$$ are all outside of the sphere and we can continue at $(a, b + 1, 0, 0, 0)$. In practice, **if outside of sphere is encountered, we set the first non-zero dimension and lower dimensions to $0$ and increase the next dimension.**

3. Note that all dimensions are symmetrical, therefore, the order of each dimension is unimportant. That is to say, for instance,
$$  (a, b, c),\ (a, c, b),\ (b, a, c),\ (b, c, a),\ (c, a, b),\ (c, b, a)
$$ have the same contribution,
$$  (a, a, b),\ (a, b, a),\ (b, a, a)
$$ have the same contribution. Then, if given a set of numbers, we can calculate all permutations at the same time. In practice, **we maintain the coordinate on each dimension is less or equal to the next dimension in our iterations, and multiply the contribution by the number of its permutations.**

## Code
```c
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

    // Number of samples for each number of dimensions
    // These numbers are tested so that the error is within 5%
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
    
    // Exact value of n-dimensional sphere
    double V_exact[101] = {0.0,
        2.0000000000,
        3.1415926536,
        4.1887902048,
        // ... Lots of data omitted
        9.4714080227e-40,
        2.3682021019e-40};
    
    mpz_t fact[N + 1];
    mpf_t V[N + 1];
    mpf_t error[N + 1];
    for (int i = 0; i <= N; i++) {
        mpf_init(V[i]);
        mpf_init(error[i]);
        mpz_init(fact[i]);
        mpz_fac_ui(fact[i], i);     // fact[i] = i!
    }
    mpf_set_ui(V[1], 2.0);
    
    mpz_t P[N + 1][N + 1];      // P[a][b] = a! / b!
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
	for (int n = 2; n <= N; n++) {      // Dimensions 2 to N
        mpz_set_ui(m, 0);               // # of points inside the sphere
        mpz_set_ui(cube, L[n] - 1);     // Volume of n-dimensional cube
        int node[n + 2];                // Coordinates
        node[n] = node[n + 1] = 0;
        for (int j = n - 1; j >= 1; j--) {
            mpz_mul_ui(cube, cube, L[n] - 1);
            node[j] = 0;
        }
        int L_sq = (L[n] - 1) * (L[n] - 1);
        while (node[n] < L[n]) {        // All points inside the cube
            int mark = 0;
            int RR = 0;
            for (int j = 1; j <= n; j++) {
                RR += node[j] * node[j];
                if (RR > L_sq) {        // Outside of the sphere
                    mark = 1;
                    break;
                }
            }
            if (mark) {                 // Outside of the sphere
                int k = 1;
                while (node[k] == 0)
                    k++;
                while (k < n && node[k] == node[k + 1])
                    k++;
                node[k]++;          // Skip unnecessary iterations
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
            int l = k;      // Obtain `l' for later use
            mpz_set(coeff, P[n][k]);        // Number of permutations
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
            mpz_mul_2exp(coeff, coeff, n - edgeCount);  // 1/2 on surfaces
            mpz_add(m, m, coeff);
            node[l]++;  // Maintain the non-decreasing order of dimensions
            memset(node, 0, l * sizeof(int));
        }
        mpf_set_z(numerator, m);
        mpf_set_z(denominator, cube);
        mpf_div(V[n], numerator, denominator);      // Volume
        
        mpf_set_d(error[n], V_exact[n]);
        mpf_div(error[n], V[n], error[n]);
        mpf_sub_ui(error[n], error[n], 1);
        mpf_abs(error[n], error[n]);
        mpf_mul_ui(error[n], error[n], 100);        // Error
        
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
```
We applied the high precision algorithm realized by the library `gmp.h`. The header `time_helper.h` is used to time the program.
```c
// time_helper.h
#ifndef TIME_HELPER_H
#define TIME_HELPER_H

#include <time.h>
#include <sys/time.h>

double cpuSecond() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

#endif
```
One of the results:
```
 1 dimension : 2.0000000000e+00
 2 dimensions: 3.0987654321e+00		Error: 1.3632%
 3 dimensions: 4.2085048011e+00		Error: 0.4707%
 4 dimensions: 4.9731748209e+00		Error: 0.7776%
 5 dimensions: 5.3049839963e+00		Error: 0.7826%
 6 dimensions: 5.2444561108e+00		Error: 1.4851%
 7 dimensions: 4.8339681901e+00		Error: 2.3113%
 8 dimensions: 4.1762454567e+00		Error: 2.8958%
 9 dimensions: 3.4028016159e+00		Error: 3.1618%
10 dimensions: 2.6337127631e+00		Error: 3.2762%
Time: 0.002293s
11 dimensions: 1.9488569443e+00		Error: 3.4368%
12 dimensions: 1.3837798311e+00		Error: 3.6335%
13 dimensions: 9.4635216533e-01		Error: 3.9229%
14 dimensions: 6.2502600332e-01		Error: 4.2988%
15 dimensions: 3.9933546614e-01		Error: 4.6907%
16 dimensions: 2.4493180345e-01		Error: 4.0799%
17 dimensions: 1.4705097002e-01		Error: 4.3054%
18 dimensions: 8.5877415063e-02		Error: 4.5426%
19 dimensions: 4.8858354273e-02		Error: 4.7977%
20 dimensions: 2.6884427860e-02		Error: 4.1754%
Time: 0.024531s
21 dimensions: 1.4561087419e-02		Error: 4.3869%
22 dimensions: 7.7097110086e-03		Error: 4.6033%
23 dimensions: 3.9943507684e-03		Error: 4.8205%
24 dimensions: 2.0110393536e-03		Error: 4.2219%
25 dimensions: 9.9984423474e-04		Error: 4.3981%
26 dimensions: 4.8763835760e-04		Error: 4.5755%
27 dimensions: 2.3346773310e-04		Error: 4.7541%
28 dimensions: 1.0980053076e-04		Error: 4.9336%
29 dimensions: 5.0386656845e-05		Error: 4.3465%
30 dimensions: 2.2901290875e-05		Error: 4.4988%
Time: 0.180790s
31 dimensions: 1.0242374617e-05		Error: 4.6514%
32 dimensions: 4.5097843107e-06		Error: 4.8039%
33 dimensions: 1.9558284356e-06		Error: 4.9564%
34 dimensions: 8.3015997791e-07		Error: 4.3957%
35 dimensions: 3.4967243218e-07		Error: 4.5268%
36 dimensions: 1.4525457711e-07		Error: 4.6582%
37 dimensions: 5.9529181050e-08		Error: 4.7896%
38 dimensions: 2.4077752462e-08		Error: 4.9211%
39 dimensions: 9.5542615740e-09		Error: 4.3927%
40 dimensions: 3.7671959437e-09		Error: 4.5070%
Time: 1.293794s
41 dimensions: 1.4673824805e-09		Error: 4.6214%
42 dimensions: 5.6480504504e-10		Error: 4.7358%
43 dimensions: 2.1488403210e-10		Error: 4.8503%
44 dimensions: 8.0830415748e-11		Error: 4.9649%
45 dimensions: 2.9890644543e-11		Error: 4.4560%
46 dimensions: 1.0997748292e-11		Error: 4.5565%
47 dimensions: 4.0036151470e-12		Error: 4.6571%
48 dimensions: 1.4423715843e-12		Error: 4.7577%
49 dimensions: 5.1436402438e-13		Error: 4.8584%
50 dimensions: 1.8160224851e-13		Error: 4.9591%
Time: 6.602065s
51 dimensions: 6.3137275164e-14		Error: 4.4741%
52 dimensions: 2.1860303258e-14		Error: 4.5631%
53 dimensions: 7.4977256965e-15		Error: 4.6522%
54 dimensions: 2.5478965457e-15		Error: 4.7413%
55 dimensions: 8.5799747580e-16		Error: 4.8305%
56 dimensions: 2.8636022763e-16		Error: 4.9197%
57 dimensions: 9.4244263964e-17		Error: 4.4604%
58 dimensions: 3.0909287152e-17		Error: 4.5398%
59 dimensions: 1.0051773179e-17		Error: 4.6192%
60 dimensions: 3.2417336704e-18		Error: 4.6987%
Time: 34.325696s
```

## Appendix

### Determine the number of samples
In order to determine the `L` for each dimension, we use a loop to test each dimension.
```c
// Compile with
// gcc -o test_size-gmp test_size-gmp.c -lgmp
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "time_helper.h"
#include <gmp.h>

#define N 200

int main(int argc, const char* argv[]) {

	if (argc < 2) {
		printf("Usage: command <number of dimensions>\n");
		printf("Usage: command <number of dimensions> <given L>\n");
		printf("Usage: command <number of dimensions> <lower bound of L> <upper bound of L>\n");
		return -1;
	}
    
	int n = 0;
	int t;
    int l = 0, r = 0;
    
    if (argc == 2) {
        t = 0;
        while (argv[1][t] != '\0') {
            n *= 10;
            n += argv[1][t] - '0';
            t++;
        }
        l = 2;
        r = 80;
    } else if (argc == 3) {
        t = 0;
        while (argv[1][t] != '\0') {
            n *= 10;
            n += argv[1][t] - '0';
            t++;
        }
        t = 0;
        while (argv[2][t] != '\0') {
            l *= 10;
            l += argv[2][t] - '0';
            t++;
        }
        r = l;
    } else if (argc == 4) {
        t = 0;
        while (argv[1][t] != '\0') {
            n *= 10;
            n += argv[1][t] - '0';
            t++;
        }
        t = 0;
        while (argv[2][t] != '\0') {
            l *= 10;
            l += argv[2][t] - '0';
            t++;
        }
        t = 0;
        while (argv[3][t] != '\0') {
            r *= 10;
            r += argv[3][t] - '0';
            t++;
        }
    }
    
    unsigned long precision = ceil(50 * log2(10) + 16);
    mpf_set_default_prec(precision);

	double V[N + 1] = {0.0,
        2.0000000000,
        3.1415926536,
        4.1887902048,
        // ...
        3.1401583155e-108,
        5.5588328420e-109};
    
    mpz_t fact[N + 1];
    for (int i = 0; i <= N; i++) {
        mpz_init(fact[i]);
        mpz_fac_ui(fact[i], i);
    }
    
    mpz_t P[N + 1][N + 1];
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= i; j++) {
            mpz_init(P[i][j]);
            mpz_tdiv_q(P[i][j], fact[i], fact[j]);
        }
    }
	
	printf("n = %d\n", n);
    
    double begin = cpuSecond();

    mpf_t Vn, error;
    mpf_init(Vn);
    mpf_init(error);
    mpz_t m;
    mpz_t cube;
    mpz_t coeff;
    mpz_init(m);
    mpz_init(cube);
    mpz_init(coeff);
    mpf_t numerator, denominator;
    mpf_init(numerator);
    mpf_init(denominator);
	int L = l - 1;
    while (++L <= r) {
        mpz_set_ui(m, 0);
        mpz_set_ui(cube, L - 1);
        int node[n + 2];
        node[n] = node[n + 1] = 0;
        for (int j = n - 1; j >= 1; j--) {
            mpz_mul_ui(cube, cube, L - 1);
            node[j] = 0;
        }
        int L_sq = (L - 1) * (L - 1);
        while (node[n] < L) {
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
                if (node[k] == 0 || node[k] == L - 1)
                    edgeCount++;
                k++;
            }
            int l = k;
            mpz_set(coeff, P[n][k]);
            int sameCount = 1;
            if (node[k] == 0 || node[k] == L - 1)
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
                    if (node[k] == 0 || node[k] == L - 1)
                        edgeCount++;
                }
                if (node[k] == 0 || node[k] == L - 1)
                    edgeCount++;
            }
            mpz_mul_2exp(coeff, coeff, n - edgeCount);
            mpz_add(m, m, coeff);
            node[l]++;
            memset(node, 0, l * sizeof(int));
        }
        mpf_set_z(numerator, m);
        mpf_set_z(denominator, cube);
        mpf_div(Vn, numerator, denominator);
        
        mpf_set_d(error, V[n]);
        mpf_div(error, Vn, error);
        mpf_sub_ui(error, error, 1);
        mpf_abs(error, error);
        mpf_mul_ui(error, error, 100);
        
        if (mpf_cmp_d(error, 5) < 0 || l == r) {
            break;
        }
    }
    printf("Time: %fs\n", cpuSecond() - begin);
    gmp_printf("V = %.*Fe\n", 10, Vn);
    printf("L = %d\n", L);
    
    if (argc == 3) {
        gmp_printf("Error: %.*Ff%%\n", 4, error);
    }

	return 0;
}
```
We used this program to calculate the volume of the 150-dimensional sphere:
```
$ ./test_size-gmp 150 21
n = 150
Time: 64.975598s
V = 8.5446136766e-73
L = 21
Error: 9.6639%
```

### Time the program
To properly time the program and obtain the average and variance, we used the following program:
```c
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
#define NUM_OF_MEASUREMENTS 20

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
        // ...
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
    
    mpz_t m;
    mpz_t cube;
    mpz_t coeff;
    mpz_init(m);
    mpz_init(cube);
    mpz_init(coeff);
    mpf_t numerator, denominator;
    mpf_init(numerator);
    mpf_init(denominator);
    
    double begin, totTime, Time[20][NUM_OF_MEASUREMENTS + 1];
    double TimeAverage[20], TimeVariance[20];
    int test = 0, num = NUM_OF_MEASUREMENTS;
    
    while (++test <= num) {
        begin = cpuSecond();
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

            if (n % 10 == 0) {
                Time[n / 10][test] = cpuSecond() - begin;
            }
        }
    }
    
    for (int i = 1; i <= N / 10; i++) {
        for (int j = 1; j <= num; j++) {
            TimeAverage[i] += Time[i][j];
            TimeVariance[i] += Time[i][j] * Time[i][j];
        }
        TimeAverage[i] /= num;
        TimeVariance[i] /= num;
        TimeVariance[i] -= TimeAverage[i] * TimeAverage[i];
        TimeVariance[i] /= (num - 1);
        TimeVariance[i] = sqrt(TimeVariance[i]);
        printf("Dimensions <= %2d: %f +- %f\n", i * 10,
               TimeAverage[i], TimeVariance[i]);
    }

	return 0;
}
```
Result:
```
Dimensions <= 10: 0.000732 +- 0.000067
Dimensions <= 20: 0.012516 +- 0.000697
Dimensions <= 30: 0.157451 +- 0.001251
Dimensions <= 40: 1.264976 +- 0.001387
Dimensions <= 50: 6.563786 +- 0.001669
Dimensions <= 60: 34.137742 +- 0.049087
```
# Exact value of n-dimensional sphere
We apply the formula:
$$  V_n = \frac{\pi^{n / 2}}{\Gamma(\frac{n}{2} + 1)}.
$$
```c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#define DIGITS 150

int main(int argc, const char* argv[]) {

	if (argc < 2) {
		printf("Usage: command <dimensions>\n");
		return -1;
	}
	int N = 0;
    int i = 0;
    while (argv[1][i] != '\0') {
        N *= 10;
        N += argv[1][i] - '0';
        i++;
    }
    if (N <= 0) {
		printf("Invalid argument. <dimensions> must be a positive integer.\n");
		return 1;
	}

	freopen("joke.txt", "w", stdout);

	unsigned long precision = ceil(DIGITS * log2(10) + 16);
	mpf_set_default_prec(precision);

	mpf_t result, pi, sqrtPi, numerator, denominator, facOverTwos, facAsT, twos, two;
	mpz_t fac;
	mpf_inits(result, sqrtPi, numerator, denominator, facOverTwos, facAsT, twos, NULL);
	mpz_init(fac);
	mpf_init_set_str(pi, "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647", 10);
	mpf_init_set_ui(two, 2);


	for (int n = 1; n <= N; n++) {
		mpf_sqrt(sqrtPi, pi);
		mpf_pow_ui(numerator, sqrtPi, n);
		if (n % 2) {		// Odd
			int m = n / 2 + 1;
			mpz_2fac_ui(fac, n);
			mpf_pow_ui(twos, two, m);
			mpf_set_z(facAsT, fac);
			mpf_div(facOverTwos, facAsT, twos);
			mpf_mul(denominator, facOverTwos, sqrtPi);
		} else {			// Even
			int m = n / 2;
			mpz_fac_ui(fac, m);
			mpf_set_z(denominator, fac);
		}
		mpf_div(result, numerator, denominator);
		if (n == N) {
            gmp_printf("%.*Fe\n", 10, result);
			break;
		}
        gmp_printf("%.*Fe,\n", 10, result);
	}

	mpf_clears(pi, result, sqrtPi, numerator, denominator, facOverTwos, facAsT, twos, two, NULL);
	mpz_clear(fac);

	return 0;
}
```
Then we have:
```c
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
```
