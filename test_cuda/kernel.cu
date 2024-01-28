
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// aerobus include
#include "../src/lib.h"
#include <stdio.h>
#include <Windows.h>

DEVICE INLINED float aexpm1(const float x) {
    using expm1 = aerobus::expm1<aerobus::i32, 9>;
    return expm1::eval<float>(x);
}

DEVICE INLINED float expm_12(const float x) {
    return aexpm1(aexpm1(aexpm1(aexpm1(
           aexpm1(aexpm1(aexpm1(aexpm1(
           aexpm1(aexpm1(aexpm1(aexpm1(x))))))))))));
}

void run_host(float* dest, const float* input, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        dest[i] = expm_12(input[i]);
    }
}

__global__ void run_device(float* dest, const float *input, size_t size)
{
    for (size_t i = threadIdx.x + blockIdx.x * blockDim.x; i < size; i += blockDim.x * gridDim.x) {
        dest[i] = expm_12(input[i]);
    }
}

float rand_range(float from, float to) {
    float s = (float) rand() / RAND_MAX;
    return from = s * (to - from);
}


int main()
{
    LARGE_INTEGER start, stop, freq;
    const size_t N = 1024*1024*128;
    float* input = (float*)malloc(N * sizeof(float));
    float* output_device = (float*)malloc(N * sizeof(float));
    float* output_host = (float*)malloc(N * sizeof(float));
    memset((void*)output_device, 0, N * sizeof(float));
    memset((void*)output_host, 0, N * sizeof(float));
    for (size_t i = 0; i < N; ++i) {
        input[i] = rand_range(-0.1F, 0.1F);
    }

    float* d_input, * d_output;
    cudaMalloc<float>(&d_input, N * sizeof(float));
    cudaMalloc<float>(&d_output, N * sizeof(float));
    cudaMemcpy((void*)d_input, input, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    ::QueryPerformanceFrequency(&freq);

    run_device<<<32 * prop.multiProcessorCount, 128>>>(d_output, d_input, N);
    cudaMemcpy((void*)output_device, (void*)d_output, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    ::QueryPerformanceCounter(&start);
    run_host(output_host, input, N);
    ::QueryPerformanceCounter(&stop);



    for (int i = 0; i < N; ++i) {
        uint32_t x_host = ((uint32_t*)output_device)[i];
        uint32_t x_device = ((uint32_t*)output_host)[i];
        if (x_host != x_device) {
            printf("error at %d : expected %.9g got %.9g", i, output_host[i], output_device[i]);
            return 1;
        }
    }

    return 0;
}
