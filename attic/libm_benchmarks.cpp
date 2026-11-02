// Benchmarks of aerobus::libm, moved out of src/benchmarks.cpp.
// aerobus::libm was removed from aerobus.h; this file does not compile against the current header.
#include "../src/aerobus.h"
#include <benchmark/benchmark.h>  // NOLINT
#include <cmath>  // NOLINT

double rand(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);  // NOLINT
}

static void BM_aero_sin_12(benchmark::State &state) {
    float *in = aerobus::aligned_malloc<float>(state.range(0), 64);
    float *out = aerobus::aligned_malloc<float>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(-0.01, 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = aerobus::libm::sin(aerobus::libm::sin(aerobus::libm::sin(
                aerobus::libm::sin(aerobus::libm::sin(aerobus::libm::sin(
                    aerobus::libm::sin(aerobus::libm::sin(aerobus::libm::sin(in[i])))))))));
        }
    }

    free(in);
    free(out);
}


static void BM_aero_fast_sin_12(benchmark::State &state) {
    float *in = aerobus::aligned_malloc<float>(state.range(0), 64);
    float *out = aerobus::aligned_malloc<float>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(-0.01, 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = aerobus::libm::fast_sin(aerobus::libm::fast_sin(aerobus::libm::fast_sin(
                aerobus::libm::fast_sin(aerobus::libm::fast_sin(aerobus::libm::fast_sin(
                    aerobus::libm::fast_sin(aerobus::libm::fast_sin(aerobus::libm::fast_sin(in[i])))))))));
        }
    }

    free(in);
    free(out);
}

static void BM_aero_cos_12(benchmark::State &state) {
    using constants = aerobus::arithmetic_helpers<float>;
    float *in = aerobus::aligned_malloc<float>(state.range(0), 64);
    float *out = aerobus::aligned_malloc<float>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(constants::pi_2() - 0.01, constants::pi_2() + 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = aerobus::libm::cos(aerobus::libm::cos(aerobus::libm::cos(
                aerobus::libm::cos(aerobus::libm::cos(aerobus::libm::cos(
                    aerobus::libm::cos(aerobus::libm::cos(aerobus::libm::cos(in[i])))))))));
        }
    }

    free(in);
    free(out);
}

static void BM_std_sin_12(benchmark::State &state) {
    double *in = aerobus::aligned_malloc<double>(state.range(0), 64);
    double *out = aerobus::aligned_malloc<double>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(-0.01, 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = ::sin(::sin(::sin(::sin(::sin(::sin(::sin(::sin(::sin(::sin(::sin(::sin(in[i]))))))))))));
        }
    }

    free(in);
    free(out);
}

static void BM_std_cos_12(benchmark::State &state) {
    using constants = aerobus::arithmetic_helpers<float>;
    double *in = aerobus::aligned_malloc<double>(state.range(0), 64);
    double *out = aerobus::aligned_malloc<double>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(constants::pi_2() - 0.01, constants::pi_2() + 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = ::cos(cos(cos(::cos(cos(cos(::cos(cos(cos(::cos(cos(cos(in[i]))))))))))));
        }
    }

    free(in);
    free(out);
}

BENCHMARK(BM_std_cos_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_aero_cos_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_std_sin_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_aero_sin_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_aero_fast_sin_12)->Range(1 << 10, 1 << 24);

BENCHMARK_MAIN();
