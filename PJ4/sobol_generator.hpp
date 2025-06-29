#ifndef SOBOL_GENERATOR_HPP
#define SOBOL_GENERATOR_HPP

#include <cstdint>

class Xorshift32 {
private:
    uint32_t state;

public:
    explicit Xorshift32(uint32_t seed = 0);
    uint32_t next();
    int32_t range(int32_t min, int32_t max);
};

void i8_sobol_generate(double* r, int m, int n, long long skip);

#endif


