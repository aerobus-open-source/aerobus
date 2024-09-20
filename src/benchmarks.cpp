#include <cstdio>
#include <cstdlib>
#include <chrono> // NOLINT
#include <cmath>
#include "./aerobus.h"

const size_t N = 134217728/8;


#ifndef DEGREE
#define DEGREE 13
#endif

template<int deg>
INLINED
double exp_12(const double x) {
    using EXP = aerobus::expm1<aerobus::i64, deg>;
    return EXP::eval(EXP::eval(EXP::eval(EXP::eval(EXP::eval(EXP::eval(
        EXP::eval(EXP::eval(EXP::eval(EXP::eval(EXP::eval(EXP::eval(x))))))))))));
}

template<size_t N, int deg>
INLINED
void vexp_12(double * in, double * out) {
    static_assert(N % 32 == 0);
    in = reinterpret_cast<double*>(__builtin_assume_aligned(in, 64));
    out = reinterpret_cast<double*>(__builtin_assume_aligned(out, 64));
#pragma omp parallel for
    for (size_t i = 0; i < N; i += 8) {
        for (size_t j = i; j < i+8; ++j) {
            out[j] = exp_12<deg>(in[j]);
        }
    }
}

INLINED
double exp_12_slow(const double x) {
    return ::expm1(::expm1(::expm1(::expm1(::expm1(::expm1(
        ::expm1(::expm1(::expm1(::expm1(::expm1(::expm1(x))))))))))));
}

template<size_t N>
INLINED
void vexp_12_slow(double * in, double * out) {
    static_assert(N % 32 == 0);
    in = reinterpret_cast<double*>(__builtin_assume_aligned(in, 64));
    out = reinterpret_cast<double*>(__builtin_assume_aligned(out, 64));
#pragma omp parallel for
    for (size_t i = 0; i < N; i += 8) {
        for (size_t j = i; j < i+8; ++j) {
            out[j] = exp_12_slow(in[j]);
        }
    }
}

double rand(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);  // NOLINT
}

template<int deg>
void run_aero(double* in, double* out_aero) {
    // warmup
    vexp_12<N, deg>(in, out_aero);
    vexp_12<N, deg>(in, out_aero);
    vexp_12<N, deg>(in, out_aero);
    double best = 1.0E9;
    for (int i = 0; i < 10; ++i) {
        auto start = std::chrono::steady_clock::now();
        vexp_12<N, deg>(in, out_aero);
        auto stop = std::chrono::steady_clock::now();
        std::chrono::duration<double> time = stop - start;
        if (time.count() < best) {
            best = time.count();
        }
    }

    printf("[aerobus deg %d] %.3e Gexp/s\n", deg, (12.0E-9 * static_cast<double>(N)) / best);
}

void verify(double* out_aero, double* out_std) {
    double err_aero = 0.0;
    double max_err = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double err = fabs(out_aero[i] - out_std[i]);
        err_aero += err;
        if (err > max_err) {
            max_err = err;
        }
    }
    err_aero /= N;
    printf("average error (vs std) : %.2e\n", err_aero);
    printf("max error (vs std) : %.2e\n", max_err);
}

int main() {
    double* in = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
    double* out_aero = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
    double* out_std = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
    double* out_vml = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
    memset(reinterpret_cast<void*>(out_aero), 0, N * sizeof(double));
    memset(reinterpret_cast<void*>(out_std), 0, N * sizeof(double));
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        in[i] = rand(-0.01, 0.01);  // pi / 6
    }

    {
        // warmup
        vexp_12_slow<N>(in, out_std);
        vexp_12_slow<N>(in, out_std);
        vexp_12_slow<N>(in, out_std);
        double best = 1.0E9;
        for (int i = 0; i < 10; ++i) {
            auto start = std::chrono::steady_clock::now();
            vexp_12_slow<N>(in, out_std);
            auto stop = std::chrono::steady_clock::now();
            std::chrono::duration<double> time = stop - start;
            if (time.count() < best) {
                best = time.count();
            }
        }

        printf("[std math] %.3e Gexp/s\n", (12.0E-9 * static_cast<double>(N)) / best);
    }

    {
        run_aero<1>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<3>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<5>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<7>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<9>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<11>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<13>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<15>(in, out_aero);
        verify(out_aero, out_std);
        run_aero<17>(in, out_aero);
        verify(out_aero, out_std);
    }

    double err_aero = 0.0;
    double max_err = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double err = fabs(out_aero[i] - out_std[i]);
        err_aero += err;
        if (err > max_err) {
            max_err = err;
        }
    }
    err_aero /= N;
    printf("average error (vs std) : %.2e\n", err_aero);
    printf("max error (vs std) : %.2e\n", max_err);
    return 0;
}
