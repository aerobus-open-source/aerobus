#include "./aerobus.h"
#include <benchmark/benchmark.h>  // NOLINT
#include <cmath>  // NOLINT

double rand(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);  // NOLINT
}

INLINED double aero_expm1_12(const double x) {
    using EXPM1 = aerobus::expm1<aerobus::i64, 13>;
    return EXPM1::eval(EXPM1::eval(EXPM1::eval(EXPM1::eval(EXPM1::eval(EXPM1::eval(
                EXPM1::eval(EXPM1::eval(EXPM1::eval(EXPM1::eval(EXPM1::eval(EXPM1::eval(x))))))))))));
}

static void BM_aero_expm1_12(benchmark::State &state) {
    double *in = aerobus::aligned_malloc<double>(state.range(0), 64);
    double *out = aerobus::aligned_malloc<double>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(-0.01, 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = aero_expm1_12(in[i]);
        }
    }

    free(in);
    free(out);
}

static void BM_std_expm1_12(benchmark::State &state) {
    float *in = aerobus::aligned_malloc<float>(state.range(0), 64);
    float *out = aerobus::aligned_malloc<float>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(-0.01, 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = ::expm1(::expm1(::expm1(::expm1(::expm1(::expm1(
                ::expm1(::expm1(::expm1(::expm1(::expm1(::expm1(in[i]))))))))))));
        }
    }

    free(in);
    free(out);
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


static void BM_aero_cos_12(benchmark::State &state) {
    using constants = aerobus::internal::arithmetic_helpers<float>;
    float *in = aerobus::aligned_malloc<float>(state.range(0), 64);
    float *out = aerobus::aligned_malloc<float>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(constants::pi_2 - 0.01, constants::pi_2 + 0.01);
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

static void BM_std_cos_12(benchmark::State &state) {
    using constants = aerobus::internal::arithmetic_helpers<float>;
    double *in = aerobus::aligned_malloc<double>(state.range(0), 64);
    double *out = aerobus::aligned_malloc<double>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(constants::pi_2 - 0.01, constants::pi_2 + 0.01);
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

static void BM_aero_hermite(benchmark::State &state) {
    double *in = aerobus::aligned_malloc<double>(state.range(0), 64);
    double *out = aerobus::aligned_malloc<double>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(-0.01, 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = aerobus::known_polynomials::hermite_phys<12>::eval<double>(in[i]);
        }
    }

    free(in);
    free(out);
}

static void BM_std_hermite(benchmark::State &state) {
    double *in = aerobus::aligned_malloc<double>(state.range(0), 64);
    double *out = aerobus::aligned_malloc<double>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(-0.01, 0.01);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = std::hermite(12, in[i]);
        }
    }

    free(in);
    free(out);
}

static void BM_compensated_horner_float(benchmark::State &state) {
    using P = aerobus::make_int_polynomial_t<aerobus::i64, 1, -11, 55, -165, 330, -462, 462, -330, 165, -55, 11, -1>;

    float *in = aerobus::aligned_malloc<float>(state.range(0), 64);
    float *out = aerobus::aligned_malloc<float>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(0.9, 1.1);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = P::compensated_eval(in[i]);
        }
    }

    free(in);
    free(out);
}
static void BM_horner_double(benchmark::State &state) {
    using P = aerobus::make_int_polynomial_t<aerobus::i64, 1, -11, 55, -165, 330, -462, 462, -330, 165, -55, 11, -1>;

    double *in = aerobus::aligned_malloc<double>(state.range(0), 64);
    double *out = aerobus::aligned_malloc<double>(state.range(0), 64);
    #pragma omp parallel for
    for (int64_t i = 0; i < state.range(0); ++i) {
        in[i] = rand(0.9, 1.1);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (int64_t i = 0; i < state.range(0); ++i) {
            out[i] = P::eval(in[i]);
        }
    }

    free(in);
    free(out);
}

BENCHMARK(BM_std_cos_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_aero_cos_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_std_sin_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_aero_sin_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_std_expm1_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_aero_expm1_12)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_std_hermite)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_aero_hermite)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_horner_double)->Range(1 << 10, 1 << 24);
BENCHMARK(BM_compensated_horner_float)->Range(1 << 10, 1 << 24);

BENCHMARK_MAIN();
