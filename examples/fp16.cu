// TO compile with nvcc -O3 -std=c++20 -arch=sm_90 fp16.cu
#include <cstdio>

#define WITH_CUDA_FP16
#include "../src/aerobus.h"

/*
change int_type to aerobus::i32 (or i64) and float_type to float (resp. double)
to see how good is the generated assembly compared to what nvcc generates for 16 bits 
*/
using int_type = aerobus::i16;
using float_type = __half2;


constexpr size_t N = 1 << 24;

template<typename T>
struct Expm1Degree;

template<>
struct Expm1Degree<float> {
    static constexpr size_t val = 11;
};

template<>
struct Expm1Degree<__half2> {
    static constexpr size_t val = 6;
};

float rand(float min, float max) {
  float range = (max - min);
  float div = RAND_MAX / range;
  return min + (rand() / div);  // NOLINT
}

template<typename T>
struct GetRandT;

template<>
struct GetRandT<float> {
    static float func(float min, float max) {
        return rand(min, max);
    }
};

template<>
struct GetRandT<__half2> {
    static __half2 func(float min, float max) {
        return __half2(__float2half(rand(min, max)), __float2half(rand(min, max)));
    }
};

using EXPM1 = aerobus::expm1<int_type, Expm1Degree<float_type>::val>;


__device__ INLINED float_type f(float_type x) {
    return EXPM1::eval(x);
}

__global__ void run(size_t N, float_type* in, float_type* out) {
    for(size_t i = threadIdx.x + blockDim.x * blockIdx.x; i < N; i += blockDim.x * gridDim.x) {
        out[i] = f(f(f(f(f(f(in[i]))))));
    }
}

int main() {
    float_type *d_in, *d_out;
    cudaMalloc<float_type>(&d_in, N * sizeof(float_type));
    cudaMalloc<float_type>(&d_out, N * sizeof(float_type));

    float_type *in = reinterpret_cast<float_type*>(malloc(N * sizeof(float_type)));
    float_type *out = reinterpret_cast<float_type*>(malloc(N * sizeof(float_type)));

    for(size_t i = 0; i < N; ++i) {
        in[i] = GetRandT<float_type>::func(-0.01, 0.01);
    }

    cudaMemcpy(d_in, in, N * sizeof(float_type), cudaMemcpyHostToDevice);

    run<<<128, 512>>>(N, d_in, d_out);

    cudaMemcpy(out, d_out, N * sizeof(float_type), cudaMemcpyDeviceToHost);

    cudaFree(d_in);
    cudaFree(d_out);
}