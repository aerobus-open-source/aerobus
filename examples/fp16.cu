// TO compile with nvcc -O3 -std=c++20 -arch=sm_90 fp16.cu
// TO GET optimal performances, modify cuda_fp16.h by adding __CUDA_FP16_CONSTEXPR__ to line 5039 (version 12.6)
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
struct Expm1Degree<double> {
    static constexpr size_t val = 18;
};

template<>
struct Expm1Degree<float> {
    static constexpr size_t val = 11;
};

template<>
struct Expm1Degree<__half2> {
    static constexpr size_t val = 6;
};

template<>
struct Expm1Degree<__half> {
    static constexpr size_t val = 6;
};

double rand(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);  // NOLINT
}

template<typename T>
struct GetRandT;

template<>
struct GetRandT<double> {
    static double func(double min, double max) {
        return rand(min, max);
    }
};

template<>
struct GetRandT<float> {
    static float func(double min, double max) {
        return (float) rand(min, max);
    }
};

template<>
struct GetRandT<__half2> {
    static __half2 func(double min, double max) {
        return __half2(__float2half((float)rand(min, max)), __float2half((float)rand(min, max)));
    }
};

template<>
struct GetRandT<__half> {
    static __half func(double min, double max) {
        return __float2half((float)rand(min, max));
    }
};

using EXPM1 = aerobus::expm1<int_type, Expm1Degree<float_type>::val>;


__device__ INLINED float_type f(float_type x) {
    return EXPM1::eval(x);
}

__global__ void run(size_t N, float_type* in, float_type* out) {
    for(size_t i = threadIdx.x + blockDim.x * blockIdx.x; i < N; i += blockDim.x * gridDim.x) {
        out[i] = f(f(f(f(f(f(f(f(f(f(f(f(in[i]))))))))))));
    }
}

#define cudaErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

int main() {
    // configure CUDA devices
    int deviceCount;
    int device = -1;
    int maxProcCount = 0;
    cudaErrorCheck(cudaGetDeviceCount(&deviceCount));
    for(int i = 0; i < deviceCount; ++i) {
        cudaDeviceProp prop;
        cudaErrorCheck(cudaGetDeviceProperties(&prop, i));
        int procCount = prop.multiProcessorCount;
        if(procCount > maxProcCount) {
            maxProcCount = procCount;
            device = i;
        }
    }

    if(device == -1) {
        ::printf("CANNOT FIND CUDA CAPABLE DEVICE -- aborting\n");
        ::abort();
    }

    cudaErrorCheck(cudaSetDevice(device));


    // allocate and populate memory
    float_type *d_in, *d_out;
    cudaErrorCheck(cudaMalloc<float_type>(&d_in, N * sizeof(float_type)));
    cudaErrorCheck(cudaMalloc<float_type>(&d_out, N * sizeof(float_type)));

    float_type *in = reinterpret_cast<float_type*>(malloc(N * sizeof(float_type)));
    float_type *out = reinterpret_cast<float_type*>(malloc(N * sizeof(float_type)));

    for(size_t i = 0; i < N; ++i) {
        in[i] = GetRandT<float_type>::func(-0.01, 0.01);
    }

    cudaErrorCheck(cudaMemcpy(d_in, in, N * sizeof(float_type), cudaMemcpyHostToDevice));

    // execute kernel and get memory back from device
    run<<<8 * maxProcCount, 256>>>(N, d_in, d_out);
    cudaErrorCheck(cudaPeekAtLastError());
    cudaErrorCheck(cudaMemcpy(out, d_out, N * sizeof(float_type), cudaMemcpyDeviceToHost));

    cudaErrorCheck(cudaFree(d_in));
    cudaErrorCheck(cudaFree(d_out));
}

// generated SASS : 

/*
HFMA2.MMA R5, R6, RZ, 0.0013885498046875, 0.0013885498046875 ;    
HFMA2 R5, R6, R5, 0.008331298828125, 0.008331298828125 ;          
HFMA2.MMA R5, R6, R5, 0.041656494140625, 0.041656494140625 ;      
HFMA2 R5, R6, R5, 0.1666259765625, 0.1666259765625 ;              
HFMA2.MMA R5, R6, R5, 0.5, 0.5 ;                                  
HFMA2 R5, R6, R5, 1, 1 ;                                          
HFMA2.MMA R5, R6, R5, RZ ;                                        
HFMA2 R7, R5, RZ.H0_H0, 0.0013885498046875, 0.0013885498046875 ;  
HFMA2.MMA R7, R5, R7, 0.008331298828125, 0.008331298828125 ;      
HFMA2 R7, R5, R7, 0.041656494140625, 0.041656494140625 ;          
HFMA2.MMA R7, R5, R7, 0.1666259765625, 0.1666259765625 ;          
HFMA2 R7, R5, R7, 0.5, 0.5 ;                                      
HFMA2.MMA R7, R5, R7, 1, 1 ;                                      
HFMA2 R7, R5, R7, RZ.H0_H0 ;                                      
HFMA2.MMA R5, R7, RZ, 0.0013885498046875, 0.0013885498046875 ;    
HFMA2 R5, R7, R5, 0.008331298828125, 0.008331298828125 ;          
HFMA2.MMA R5, R7, R5, 0.041656494140625, 0.041656494140625 ;      
HFMA2 R5, R7, R5, 0.1666259765625, 0.1666259765625 ;              
HFMA2.MMA R5, R7, R5, 0.5, 0.5 ;                                  
HFMA2 R5, R7, R5, 1, 1 ;                                          
HFMA2.MMA R5, R7, R5, RZ ;                                        
HFMA2 R6, R5, RZ.H0_H0, 0.0013885498046875, 0.0013885498046875 ;  
HFMA2.MMA R6, R5, R6, 0.008331298828125, 0.008331298828125 ;      
HFMA2 R6, R5, R6, 0.041656494140625, 0.041656494140625 ;          
HFMA2.MMA R6, R5, R6, 0.1666259765625, 0.1666259765625 ;          
HFMA2 R6, R5, R6, 0.5, 0.5 ;                                      
HFMA2.MMA R6, R5, R6, 1, 1 ;                                      
HFMA2 R6, R5, R6, RZ.H0_H0 ;                                      
HFMA2.MMA R5, R6, RZ, 0.0013885498046875, 0.0013885498046875 ;    
HFMA2 R5, R6, R5, 0.008331298828125, 0.008331298828125 ;          
HFMA2.MMA R5, R6, R5, 0.041656494140625, 0.041656494140625 ;      
HFMA2 R5, R6, R5, 0.1666259765625, 0.1666259765625 ;              
HFMA2.MMA R5, R6, R5, 0.5, 0.5 ;                                  
HFMA2 R5, R6, R5, 1, 1 ;                                          
HFMA2.MMA R5, R6, R5, RZ ;                                        
HFMA2 R6, R5, RZ.H0_H0, 0.0013885498046875, 0.0013885498046875 ;  
HFMA2.MMA R6, R5, R6, 0.008331298828125, 0.008331298828125 ;      
HFMA2 R6, R5, R6, 0.041656494140625, 0.041656494140625 ;          
HFMA2.MMA R6, R5, R6, 0.1666259765625, 0.1666259765625 ;          
HFMA2 R6, R5, R6, 0.5, 0.5 ;                                      
HFMA2.MMA R6, R5, R6, 1, 1 ;                                      
HFMA2 R6, R5, R6, RZ.H0_H0 ;                                      
HFMA2.MMA R5, R6, RZ, 0.0013885498046875, 0.0013885498046875 ;    
HFMA2 R5, R6, R5, 0.008331298828125, 0.008331298828125 ;          
HFMA2.MMA R5, R6, R5, 0.041656494140625, 0.041656494140625 ;      
HFMA2 R5, R6, R5, 0.1666259765625, 0.1666259765625 ;              
HFMA2.MMA R5, R6, R5, 0.5, 0.5 ;                                  
HFMA2 R5, R6, R5, 1, 1 ;                                          
HFMA2.MMA R5, R6, R5, RZ ;                                        
HFMA2 R6, R5, RZ.H0_H0, 0.0013885498046875, 0.0013885498046875 ;  
HFMA2.MMA R6, R5, R6, 0.008331298828125, 0.008331298828125 ;      
HFMA2 R6, R5, R6, 0.041656494140625, 0.041656494140625 ;          
HFMA2.MMA R6, R5, R6, 0.1666259765625, 0.1666259765625 ;          
HFMA2 R6, R5, R6, 0.5, 0.5 ;                                      
HFMA2.MMA R6, R5, R6, 1, 1 ;                                      
HFMA2 R6, R5, R6, RZ.H0_H0 ;                                      
HFMA2.MMA R5, R6, RZ, 0.0013885498046875, 0.0013885498046875 ;    
HFMA2 R5, R6, R5, 0.008331298828125, 0.008331298828125 ;          
HFMA2.MMA R5, R6, R5, 0.041656494140625, 0.041656494140625 ;      
HFMA2 R5, R6, R5, 0.1666259765625, 0.1666259765625 ;              
HFMA2.MMA R5, R6, R5, 0.5, 0.5 ;                                  
HFMA2 R5, R6, R5, 1, 1 ;                                          
HFMA2.MMA R5, R6, R5, RZ ;                                        
HFMA2 R6, R5, RZ.H0_H0, 0.0013885498046875, 0.0013885498046875 ;  
HFMA2.MMA R6, R5, R6, 0.008331298828125, 0.008331298828125 ;      
HFMA2 R6, R5, R6, 0.041656494140625, 0.041656494140625 ;          
HFMA2.MMA R6, R5, R6, 0.1666259765625, 0.1666259765625 ;          
HFMA2 R6, R5, R6, 0.5, 0.5 ;                                      
HFMA2.MMA R6, R5, R6, 1, 1 ;                                      
HFMA2 R6, R5, R6, RZ.H0_H0 ;                                      
HFMA2.MMA R5, R6, RZ, 0.0013885498046875, 0.0013885498046875 ;    
HFMA2 R5, R6, R5, 0.008331298828125, 0.008331298828125 ;          
HFMA2.MMA R5, R6, R5, 0.041656494140625, 0.041656494140625 ;      
HFMA2 R5, R6, R5, 0.1666259765625, 0.1666259765625 ;              
HFMA2.MMA R5, R6, R5, 0.5, 0.5 ;                                  
HFMA2 R5, R6, R5, 1, 1 ;                                          
HFMA2.MMA R6, R6, R5, RZ ;                                        
HFMA2 R5, R6, RZ.H0_H0, 0.0013885498046875, 0.0013885498046875 ;  
HFMA2.MMA R5, R6, R5, 0.008331298828125, 0.008331298828125 ;      
HFMA2 R5, R6, R5, 0.041656494140625, 0.041656494140625 ;          
HFMA2.MMA R5, R6, R5, 0.1666259765625, 0.1666259765625 ;          
HFMA2 R5, R6, R5, 0.5, 0.5 ;                                      
HFMA2.MMA R7, R6, R5, 1, 1 ;                                      
IADD3.X R5, R8, UR11, RZ, P0, !PT ;                               
IADD3 R3, P0, R2, R3, RZ ;                                        
IADD3.X R0, RZ, R0, RZ, P0, !PT ;                                 
ISETP.GE.U32.AND P0, PT, R3, UR8, PT ;                            
HFMA2 R7, R6, R7, RZ.H0_H0 ;                                      
ISETP.GE.U32.AND.EX P0, PT, R0, UR9, PT, P0 ;                     
STG.E desc[UR6][R4.64], R7 ;                  
*/