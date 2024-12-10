#include <hip_runtime.h>
#include "./utilities.h"
#include "../src/aerobus.h"

// to compile with hipcc -std=c++20 -Wno-unused-result -I /usr/include/hip/ -O3 libm_hip.cpp

__global__ void runlibm(float *out, float *in, int N) {
    for (int i = threadIdx.x + blockDim.x * blockIdx.x; i < N; i += blockDim.x * gridDim.x) {
        out[i] = sin(in[i]);
    }
}

__global__ void runaero(float* out, float *in, int N) {
    for (int i = threadIdx.x + blockDim.x * blockIdx.x; i < N; i += blockDim.x * gridDim.x) {
        out[i] = aerobus::libm::sin(in[i]);
    }
}

static constexpr size_t N = 1 << 20;

int main() {
    hipEvent_t start, stop;
    hipEventCreate(&start);
    hipEventCreate(&stop);
    float milliseconds = 0;
    int deviceCount;
    int device = -1;
    int maxProcCount = 0;
    hipErrorCheck(hipGetDeviceCount(&deviceCount));
    for (int i = 0; i < deviceCount; ++i) {
        hipDeviceProp_t prop;
        hipErrorCheck(hipGetDeviceProperties(&prop, i));
        int procCount = prop.multiProcessorCount;
        if (procCount > maxProcCount) {
            maxProcCount = procCount;
            device = i;
        }
    }

    if (device == -1) {
        ::printf("CANNOT FIND HIP CAPABLE DEVICE -- aborting\n");
        ::abort();
    }

    hipErrorCheck(hipSetDevice(device));

    // memory allocation
    float *in, *d_in, *out, *d_out;
    float *rout, *rd_out;
    hipErrorCheck(hipMalloc<float>(&d_in, N * sizeof(float)));
    hipErrorCheck(hipMalloc<float>(&d_out, N * sizeof(float)));
    hipErrorCheck(hipMalloc<float>(&rd_out, N * sizeof(float)));
    in = reinterpret_cast<float*>(malloc(N * sizeof(float)));
    out = reinterpret_cast<float*>(malloc(N * sizeof(float)));
    rout = reinterpret_cast<float*>(malloc(N * sizeof(float)));

    // populate
    for (size_t i = 0; i < N; ++i) {
        in[i] = GetRandT<float>::func(-10, 10);
    }

    hipErrorCheck(hipMemcpy(d_in, in, N * sizeof(float), hipMemcpyHostToDevice));

    // execute kernel and get memory back from device

    hipEventRecord(start);
    runlibm<<<128, 256>>>(rd_out, d_in, N);
    hipErrorCheck(hipPeekAtLastError());
    hipEventRecord(stop);
    hipEventSynchronize(stop);
    hipEventElapsedTime(&milliseconds, start, stop);
    ::printf("reference time : %f\n", milliseconds);


    hipEventRecord(start);
    runaero<<<128, 256>>>(d_out, d_in, N);
    hipErrorCheck(hipPeekAtLastError());
    hipEventRecord(stop);
    hipEventSynchronize(stop);
    hipEventElapsedTime(&milliseconds, start, stop);
    ::printf("aero time : %f\n", milliseconds);

    hipErrorCheck(hipMemcpy(out, d_out, N * sizeof(float), hipMemcpyDeviceToHost));
    hipErrorCheck(hipMemcpy(rout, rd_out, N * sizeof(float), hipMemcpyDeviceToHost));

    hipErrorCheck(hipFree(d_in));
    hipErrorCheck(hipFree(d_out));

    // verify results
    int erro_count = 0;
    int ok_count = 0;
    float error = 0;
    FILE* f = fopen("errors.txt", "w+");
    for (int i = 0; i < N; ++i) {
        if (out[i] > rout[i]) {
            erro_count += 1;
            error += out[i] - rout[i];
            fprintf(f, "too large -- input = %.12f -- expected %.12f got %.12f\n",
                in[i], rout[i], out[i]);
        } else if (out[i] < rout[i]) {
            erro_count += 1;
            error += rout[i] - out[i];
            fprintf(f, "too small -- input = %.12f -- expected %.12f got %.12f\n",
                in[i], rout[i], out[i]);
        } else {
            ok_count += 1;
        }
    }
    fclose(f);

    ::printf("errors : %d, OK : %d\n", erro_count, ok_count);
    ::printf("average error : %f\n", error / N);

    return 0;
}