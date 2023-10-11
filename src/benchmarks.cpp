#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include "lib.h"
#include "vml_math.h"

const size_t N = 134217728 / 4;

INLINED
double sin_12(const double x) {
	using V = aerobus::sin<aerobus::i64, 13>;
	return V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(
		V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(x))))))))))));
}

template<size_t N>
INLINED
void vsin_12(double * in, double * out) {
	static_assert(N % 32 == 0);
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
#pragma omp parallel for
	for (int i = 0; i < N; i += 8) {
		#pragma simd
		for(int j = i; j < i+8; ++j) {
			out[j] = sin_12(in[j]);
		}
	}
}

INLINED
double sin_12_slow(const double x) {
	return ::sin(::sin(::sin(::sin(::sin(::sin(
		::sin(::sin(::sin(::sin(::sin(::sin(x))))))))))));
}

template<size_t N>
INLINED
void vsin_slow(double * in, double * out) {
	static_assert(N % 32 == 0);
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
#pragma omp parallel for
	for (int i = 0; i < N; i += 8) {
		#pragma simd
		for(int j = i; j < i+8; ++j) {
			out[j] = sin_12_slow(in[j]);
		}
	}
}


template<size_t N>
INLINED
void vml_vsin_12(double * in, double * out) {
	static_assert(N % 32 == 0);
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
	vml_vsin(in, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
	vml_vsin(out, out, (unsigned) N);
}

int main() {
	double* in = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	double* out = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	for(int i = 0; i < N; ++i) {
		in[i] = 0.1;
	}

	{
		// warmup
		vsin_12<N>(in, out);
		vsin_12<N>(in, out);
		vsin_12<N>(in, out);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			vsin_12<N>(in, out);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[aerobus] time for %zu sin  (compound 12 times) : %lf\n", N, best);
		printf("[aerobus] achieved Gflops : %lf\n", ((double)N * 13.0 * 12.0E-9 / best));
	}

	{
		// warmup
		vsin_slow<N>(in, out);
		vsin_slow<N>(in, out);
		vsin_slow<N>(in, out);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			vsin_slow<N>(in, out);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[std::math] time for %zu sin (compound 12 times) : %lf\n", N, best);
		printf("[std::math] achieved Gflops : %lf\n", ((double)N * 13.0 * 12.0E-9 / best));
	}
	{
		// warmup
		vml_vsin_12<N>(in, out);
		vml_vsin_12<N>(in, out);
		vml_vsin_12<N>(in, out);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			vml_vsin_12<N>(in, out);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[vml] time for %zu sin  (compound 12 times) : %lf\n", N, best);
		printf("[vml] achieved Gflops : %lf\n", ((double)N * 13.0 * 12.0E-9 / best));
	}
}