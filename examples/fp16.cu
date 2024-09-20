// TO compile with nvcc -O3 -std=c++20 -arch=sm_90 fp16.cu
#include <cstdio>

#define WITH_CUDA_FP16
#include "../src/aerobus.h"

using P = aerobus::known_polynomials::hermite_phys<8, aerobus::i16>;
using P2 = aerobus::known_polynomials::hermite_phys<8, aerobus::i64>;


__device__ __forceinline__ __half f(__half x) {
    return P::eval(x);
}

__global__ void run(size_t N, __half* in, __half* out) {
    for(size_t i = threadIdx.x + blockDim.x * blockIdx.x; i < N; i += blockDim.x * gridDim.x) {
        out[i] = f(in[i]);
    }
}

int main() {
    printf("hermite<8, int16>(1) == %f\n", __half2float(P::eval<__half>(__float2half(1.0F))));
    printf("hermite<8, int64>(1) == %f\n", P2::eval(1.0));
}