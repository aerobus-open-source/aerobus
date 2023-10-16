#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include "lib.h"
#include <immintrin.h>

const size_t N = 134217728 / 16;


extern "C" __m512d _ZGVeN8v_sin(__m512d x);

#ifndef DEGREE 
#define DEGREE 13
#endif

using AERO_SIN = aerobus::sin<aerobus::i64, DEGREE>;
INLINED
double sin_12(const double x) {

	return AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(
		AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(x))))))))))));
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
void fast_math_vsin_12(double * in, double * out) {
	static_assert(N % 32 == 0);
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
	#pragma omp parallel for
	for (unsigned int i = 0; i < N; i += 8)
	{
		_mm512_store_pd((void*)(out + i), 
			_ZGVeN8v_sin(_ZGVeN8v_sin(_ZGVeN8v_sin(_ZGVeN8v_sin(
			_ZGVeN8v_sin(_ZGVeN8v_sin(_ZGVeN8v_sin(_ZGVeN8v_sin(
			_ZGVeN8v_sin(_ZGVeN8v_sin(_ZGVeN8v_sin(_ZGVeN8v_sin(
				_mm512_load_pd((void*)(in + i)))))))))))))));
	}
}

double rand(double min, double max)
{
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

int main() {
	printf("RUN WITH TAYLOR EXPANSION OF DEGREE %d\n", DEGREE);
	printf("USED POLYNOMIAL AS SIN APPROX : %s\n", AERO_SIN::to_string().c_str());
	double* in = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	double* out_aero = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	double* out_std = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	double* out_vml = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	memset((void*) out_aero, 0, N * sizeof(double));
	memset((void*) out_std, 0, N * sizeof(double));
	memset((void*) out_vml, 0, N * sizeof(double));
	#pragma omp parallel for
	for(int i = 0; i < N; ++i) {
		in[i] = rand(-0.5, 0.5); // pi / 6
	}

	{
		// warmup
		vsin_12<N>(in, out_aero);
		vsin_12<N>(in, out_aero);
		vsin_12<N>(in, out_aero);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			vsin_12<N>(in, out_aero);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[aerobus] %.3e Gsin/s\n", 12.0E-9 * (double) N / best);
	}

	{
		// warmup
		vsin_slow<N>(in, out_std);
		vsin_slow<N>(in, out_std);
		vsin_slow<N>(in, out_std);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			vsin_slow<N>(in, out_std);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[std math] %.3e Gsin/s\n", 12.0E-9 * (double) N / best);
	}
	{
		// warmup
		fast_math_vsin_12<N>(in, out_vml);
		fast_math_vsin_12<N>(in, out_vml);
		fast_math_vsin_12<N>(in, out_vml);
		double best = 1.0E9;
		for (int i = 0; i < 10; ++i) {
			auto start = std::chrono::steady_clock::now();
			fast_math_vsin_12<N>(in, out_vml);
			auto stop = std::chrono::steady_clock::now();
			std::chrono::duration<double> time = stop - start;
			if (time.count() < best) {
				best = time.count();
			}
		}

		printf("[std fast math] %.3e Gsin/s\n", 12.0E-9 * (double) N / best);
	}

	printf("VERIFY :\n");
	double err_aero = 0.0;
	double max_err = 0.0;
	for(int i = 0; i < N; ++i) {
		double err = fabs(out_aero[i] - out_std[i]);
		err_aero += err;
		if (err > max_err) {
			max_err = err;
		}
	}
	err_aero /= (double) N;
	printf("average error for aerobus (compared to std) : %.2e\n", err_aero);
	printf("max error for aerobus (compared to std) : %.2e\n", max_err);
	printf("OK\n");
	return 0;
}