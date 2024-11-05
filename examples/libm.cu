#define WITH_CUDA_FP16
#define WITH_LIBM
#include "../src/aerobus.h"
#include "./utilities.h"

// must be compiled using --expt-relaxed-constexpr or it generates crap (don't know why)

__global__ void runfloat(__half *reference_d, __half *reference_u, __half *in, int N) {
    for(int i = threadIdx.x + blockDim.x * blockIdx.x; i < N; i += blockDim.x * gridDim.x) {
        float x = sin(__half2float(in[i]));
        reference_u[i] = __float2half_ru(x);
        reference_d[i] = __float2half_rd(x);
    }
}

__global__ void runhalf(__half* out, __half *in, int N) {
    for(int i = threadIdx.x + blockDim.x * blockIdx.x; i < N; i += blockDim.x * gridDim.x) {
        out[i] = aerobus::libm::sin(in[i]);
    }
}

static constexpr size_t N = 1 << 22;

int main() {  // configure CUDA devicescudaEvent_t start, stop;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float milliseconds = 0;
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
    int blockSizef;   // The launch configurator returned block size 
    int minGridSizef; // The minimum grid size needed to achieve the 
                    // maximum occupancy for a full device launch 

    cudaErrorCheck(cudaOccupancyMaxPotentialBlockSize( &minGridSizef, &blockSizef, &runfloat, 0, 0)); 

    // memory allocation
    __half *in, *d_in, *out, *d_out;
    __half *reference_u, *d_reference_u;
    __half *reference_d, *d_reference_d;
    cudaErrorCheck(cudaMalloc<__half>(&d_in, N * sizeof(__half)));
    cudaErrorCheck(cudaMalloc<__half>(&d_out, N * sizeof(__half)));
    cudaErrorCheck(cudaMalloc<__half>(&d_reference_d, N * sizeof(__half)));
    cudaErrorCheck(cudaMalloc<__half>(&d_reference_u, N * sizeof(__half)));
    in = reinterpret_cast<__half*>(malloc(N * sizeof(__half)));
    out = reinterpret_cast<__half*>(malloc(N * sizeof(__half)));
    reference_u = reinterpret_cast<__half*>(malloc(N * sizeof(__half)));
    reference_d = reinterpret_cast<__half*>(malloc(N * sizeof(__half)));

    // populate
    for(size_t i = 0; i < N; ++i) {
        in[i] = GetRandT<__half>::func(-10, 10);
    }

    cudaErrorCheck(cudaMemcpy(d_in, in, N * sizeof(__half), cudaMemcpyHostToDevice));
  
    // execute kernel and get memory back from device
    cudaEventRecord(start);
    runhalf<<<minGridSizef, blockSizef>>>(d_out, d_in, N);
    cudaErrorCheck(cudaPeekAtLastError());
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&milliseconds, start, stop);
    ::printf("reference time : %f\n", milliseconds);


    cudaEventRecord(start);
    runfloat<<<minGridSizef, blockSizef>>>(d_reference_d, d_reference_u, d_in, N);
    cudaErrorCheck(cudaPeekAtLastError());
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&milliseconds, start, stop);
    ::printf("aero time : %f\n", milliseconds);
    cudaErrorCheck(cudaMemcpy(out, d_out, N * sizeof(__half), cudaMemcpyDeviceToHost));
    cudaErrorCheck(cudaMemcpy(reference_u, d_reference_u, N * sizeof(__half), cudaMemcpyDeviceToHost));
    cudaErrorCheck(cudaMemcpy(reference_d, d_reference_d, N * sizeof(__half), cudaMemcpyDeviceToHost));

    cudaErrorCheck(cudaFree(d_in));
    cudaErrorCheck(cudaFree(d_out));
    cudaErrorCheck(cudaFree(d_reference_d));
    cudaErrorCheck(cudaFree(d_reference_u));

    // verify results
    int erro_count = 0;
    int ok_count = 0;
    float error = 0;
    FILE* f = fopen("errors.txt", "w+");
    for(int i = 0; i < N; ++i) {
        if(out[i] > reference_u[i]) {
            erro_count += 1;
            error += __half2float(out[i] - reference_u[i]);
            fprintf(f, "too large -- input = %.12f -- expected %.12f to be in [%.12f, %.12f]\n", 
                __half2float(in[i]), __half2float(out[i]), 
                __half2float(reference_d[i]), __half2float(reference_u[i]));
        } else if(out[i] < reference_d[i]) {
            erro_count += 1;
            error += __half2float(reference_d[i] - out[i]);
            fprintf(f, "too small -- input = %.12f -- expected %.12f to be in [%.12f, %.12f]\n", 
                __half2float(in[i]), __half2float(out[i]), 
                __half2float(reference_d[i]), __half2float(reference_u[i]));
        } else {
            ok_count += 1;
        }   
    }
    fclose(f);

    ::printf("errors : %d, OK : %d\n", erro_count, ok_count);
    ::printf("average error : %f\n", error / N);

    return 0;
}
