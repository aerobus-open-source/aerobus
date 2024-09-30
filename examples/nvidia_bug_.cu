// compile with nvcc -O3 -std=c++20 -arch=sm_90 --use_fast_math --keep 
// to see how bad is the assembly for fp16 compared to the fp32 version

#include <cstdint>
#include <cuda_fp16.h>

template<typename Out, typename In>
struct cast_op {
    template<In x>
    static constexpr __host__ __device__ __forceinline__ Out f() {
        return static_cast<Out>(x);
    }
};

__global__ void run__half(__half* ptr) {
    *ptr = cast_op<__half, int16_t>::template f<2>() / cast_op<__half, int16_t>::template f<3>();
}

__global__ void run__float(float* ptr) {
    *ptr = cast_op<float, int32_t>::template f<2>() / cast_op<float, int32_t>::template f<3>();
}

int main() {

}
