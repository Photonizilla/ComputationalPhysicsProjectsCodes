// Compile to library:
// nvcc -Xcompiler -fPIC -c libFractalGPU.cu -o libFractalGPU.o -lm
// nvcc -shared -Xcompiler -fPIC libFractalGPU.o -o libFractalGPU.so -lcudart
//

#include <math.h>
#include <cuComplex.h>
#include "error_check.h"
#include <stdio.h>

__device__ static const unsigned char listR[] = {0, 137, 76, 137, 20, 66, 88, 116, 115, 168, 203, 255, 248, 240, 238, 235, 136};
__device__ static const unsigned char listG[] = {0, 43, 31, 61, 50, 148, 191, 249, 246, 247, 250, 251, 214, 153, 110, 64, 33};
__device__ static const unsigned char listB[] = {0, 142, 141, 246, 245, 247, 249, 253, 156, 77, 80, 84, 72, 55, 43, 37, 17};

__global__ static void rootAtPixel(unsigned char* rootResult, const int N, const double Lx, const double Rx, const double Ly, const double Ry) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    if (col >= N || row >= N)
        return;
    cuDoubleComplex Z = make_cuDoubleComplex((Rx - Lx) / (N - 1) * col + Lx, (Ry - Ly) / (N - 1) * row + Ly);
    cuDoubleComplex solution1 = make_cuDoubleComplex(1, 0);
    cuDoubleComplex solution2 = make_cuDoubleComplex(-0.5, sqrt(0.75));
    cuDoubleComplex solution3 = make_cuDoubleComplex(-0.5, -sqrt(0.75));
    const double r = (Rx - Lx < Ry - Ly) ? (Rx - Lx) : (Ry - Ly);
    while (cuCabs(cuCsub(Z, solution1)) > r / N && cuCabs(cuCsub(Z, solution2)) > r / N && cuCabs(cuCsub(Z, solution3)) > r / N) {
        Z = cuCsub(Z, cuCdiv(cuCsub(cuCmul(cuCmul(Z, Z), Z), make_cuDoubleComplex(1, 0)), cuCmul(cuCmul(make_cuDoubleComplex(3, 0), Z), Z)));
    }
    double imagPart = cuCimag(Z);
    if (imagPart < -0.5) {
        rootResult[(row * N + col) * 3 + 0] = 85;
        rootResult[(row * N + col) * 3 + 1] = 177;
        rootResult[(row * N + col) * 3 + 2] = 71;
    } else if (imagPart > 0.5) {
        rootResult[(row * N + col) * 3 + 0] = 63;
        rootResult[(row * N + col) * 3 + 1] = 132;
        rootResult[(row * N + col) * 3 + 2] = 247;
    } else {
        rootResult[(row * N + col) * 3 + 0] = 193;
        rootResult[(row * N + col) * 3 + 1] = 76;
        rootResult[(row * N + col) * 3 + 2] = 66;
    }
}

__global__ static void orderAtPixel(unsigned char* orderResult, const int N, const double Lx, const double Rx, const double Ly, const double Ry) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    if (col >= N || row >= N)
        return;
    cuDoubleComplex Z = make_cuDoubleComplex((Rx - Lx) / (N - 1) * col + Lx, (Ry - Ly) / (N - 1) * row + Ly);
    unsigned char count = 0;
    cuDoubleComplex solution1 = make_cuDoubleComplex(1, 0);
    cuDoubleComplex solution2 = make_cuDoubleComplex(-0.5, sqrt(0.75));
    cuDoubleComplex solution3 = make_cuDoubleComplex(-0.5, -sqrt(0.75));
    const double r = (Rx - Lx < Ry - Ly) ? (Rx - Lx) : (Ry - Ly);
    while (cuCabs(cuCsub(Z, solution1)) > r / N && cuCabs(cuCsub(Z, solution2)) > r / N && cuCabs(cuCsub(Z, solution3)) > r / N) {
        Z = cuCsub(Z, cuCdiv(cuCsub(cuCmul(cuCmul(Z, Z), Z), make_cuDoubleComplex(1, 0)), cuCmul(cuCmul(make_cuDoubleComplex(3, 0), Z), Z)));
        if (count < 255)
            count++;
    }
    if (count <= 16) {
        orderResult[(row * N + col) * 3 + 0] = listR[count];
        orderResult[(row * N + col) * 3 + 1] = listG[count];
        orderResult[(row * N + col) * 3 + 2] = listB[count];
    } else {
        orderResult[(row * N + col) * 3 + 0] = (unsigned char)(136.0 * exp(-0.0625 * (count - 16)));
        orderResult[(row * N + col) * 3 + 1] = (unsigned char)(33.0 * exp(-0.0625 * (count - 16)));
        orderResult[(row * N + col) * 3 + 2] = (unsigned char)(17.0 * exp(-0.0625 * (count - 16)));
    }
}

extern "C" {

void fractalRootCalc(unsigned char* rootOutput, const int N, const double Lx, const double Rx, const double Ly, const double Ry) {
    unsigned char* rootColorsDevice = NULL;
    cudaMalloc((void**)&rootColorsDevice, N * N * 3 * sizeof(unsigned char));
    CHECK(cudaGetLastError());
    
    dim3 dimGrid(ceil(N / 16.0), ceil(N / 16.0), 1);
    dim3 dimBlock(16, 16, 1);
    printf("GPU Started...\n");
    rootAtPixel<<<dimGrid, dimBlock>>>(rootColorsDevice, N, Lx, Rx, Ly, Ry);
    CHECK(cudaDeviceSynchronize());
    CHECK(cudaGetLastError());
    
    cudaMemcpy(rootOutput, rootColorsDevice, N * N * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    CHECK(cudaGetLastError());
    
    cudaFree(rootColorsDevice);
    rootColorsDevice = NULL;
    CHECK(cudaGetLastError());
    CHECK(cudaDeviceReset());
    printf("fractalCalc Completed!\n");
}

void fractalOrderCalc(unsigned char* orderOutput, const int N, const double Lx, const double Rx, const double Ly, const double Ry) {
    unsigned char* orderColorsDevice = NULL;
    cudaMalloc((void**)&orderColorsDevice, N * N * 3 * sizeof(unsigned char));
    CHECK(cudaGetLastError());
    
    dim3 dimGrid(ceil(N / 16.0), ceil(N / 16.0), 1);
    dim3 dimBlock(16, 16, 1);
    printf("GPU Started...\n");
    orderAtPixel<<<dimGrid, dimBlock>>>(orderColorsDevice, N, Lx, Rx, Ly, Ry);
    CHECK(cudaDeviceSynchronize());
    CHECK(cudaGetLastError());
    
    cudaMemcpy(orderOutput, orderColorsDevice, N * N * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    CHECK(cudaGetLastError());
    
    cudaFree(orderColorsDevice);
    orderColorsDevice = NULL;
    CHECK(cudaGetLastError());
    CHECK(cudaDeviceReset());
    printf("fractalCalc Completed!\n");
}

}

