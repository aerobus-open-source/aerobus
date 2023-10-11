#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include "lib.h"

const size_t N = 134217728 / 4;

INLINED
double expm1_12(const double x) {
	using V = aerobus::expm1<aerobus::i64, 13>;
	return V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(
		V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(x))))))))))));
}

template<size_t N>
INLINED
void vexpm1_12(double * in, double * out) {
	static_assert(N % 32 == 0);
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
#pragma omp parallel for
	for (int i = 0; i < N; i += 8) {
		#pragma simd
		for(int j = i; j < i+8; ++j) {
			out[j] = expm1_12(in[j]);
		}
	}
}

INLINED
double expm1_12_slow(const double x) {
	return ::expm1(::expm1(::expm1(::expm1(::expm1(::expm1(
		::expm1(::expm1(::expm1(::expm1(::expm1(::expm1(x))))))))))));
}

template<size_t N>
INLINED
void vexpm1_slow(double * in, double * out) {
	static_assert(N % 32 == 0);
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
#pragma omp parallel for
	for (int i = 0; i < N; i += 8) {
		#pragma simd
		for(int j = i; j < i+8; ++j) {
			out[j] = expm1_12_slow(in[j]);
		}
	}
}

int main() {
	double* in = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	double* out = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	for(int i = 0; i < N; ++i) {
		in[i] = 0.1;
	}

	{
		// warmup
		vexpm1_12<N>(in, out);
		vexpm1_12<N>(in, out);
		vexpm1_12<N>(in, out);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			vexpm1_12<N>(in, out);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[aerobus] time for %zu exp12 : %lf\n", N, best);
		printf("[aerobus] achieved Gflops : %lf\n", ((double)N * 13.0 * 12.0E-9 / best));
	}

	{
		// warmup
		vexpm1_slow<N>(in, out);
		vexpm1_slow<N>(in, out);
		vexpm1_slow<N>(in, out);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			vexpm1_slow<N>(in, out);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[std::math] time for %zu exp12 : %lf\n", N, best);
		printf("[std::math] achieved Gflops : %lf\n", ((double)N * 13.0 * 12.0E-9 / best));
	}
}