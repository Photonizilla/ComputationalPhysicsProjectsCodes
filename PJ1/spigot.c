#include <stdio.h>
#include <math.h>
#include <time.h>
#define MAX_C 3500
int a=10000, b, c=MAX_C, d, e, f[MAX_C+1], g;
int main(int argc, char* argv[]) {
    clock_t start = clock();
    for(;b-c;)
        f[b++]=a/5;
    for(; d=0, g=c*2; c -= 14, printf("%.4d",e+d/a), e=d%a)
        for(b=c; d+=f[b]*a, f[b]=d%--g, d/=g--, --b; d*=b);
    clock_t stop = clock();
    printf("\nTime: %fs\n", (stop - start) / (float)CLOCKS_PER_SEC);
    return 0;
}
