#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include "lib.h"
#include <immintrin.h>

constexpr int64_t N = 32*1024*1024;


// sometimes, gcc and clang don't use fast math lib calls and stay scalar
// using that, we force them to use vectorized version, no matter the optimization flags
// unfortunately, gcc doesn't implement _mmXXX_sin_pd so we have to rely on extern libm functions
#ifdef __AVX512F__
#ifdef __MSC_VER
extern "C" __m512d __vdecl_sin8(__m512d x);
INLINED __m512d vec_sin(__m512d x) {
	return __vdecl_sin8(x);
}
#else
extern "C" __m512d _ZGVeN8v_sin(__m512d x); 
INLINED __m512d vec_sin(__m512d x) {
	return _ZGVeN8v_sin(x);
}
#endif
#else
	#ifdef __AVX2__
#ifdef _MSC_VER
extern "C" __m256d __vdecl_sin4(__m256d x);
INLINED __m256d vec_sin(__m256d x) {
	return __vdecl_sin4(x);
}
#else 
extern "C" __m256d _ZGVdN4v_sin(__m256d x);
INLINED __m256d vec_sin(__m256d x) {
	return _ZGVdN4v_sin(x);
}
#endif
	#else
#error "benchmark require avx2 or avx512"
	#endif
#endif



template<int deg>
INLINED double std_laguerre(const double x) {
	return 1.0 - std::laguerre(deg, x);
}

template<int deg>
INLINED double aero_laguerre(const double x) {
	// 1 - L_12
	using H = typename aerobus::pq64::template sub_t<
		typename aerobus::pq64::one,
		typename aerobus::known_polynomials::laguerre<deg>>;
	return H::template eval<double>(x);
}

template<int deg>
INLINED double std_laguerre_12(const double x) {
	return std_laguerre<deg>(std_laguerre<deg>(std_laguerre<deg>(std_laguerre<deg>(
		std_laguerre<deg>(std_laguerre<deg>(std_laguerre<deg>(std_laguerre<deg>(
			std_laguerre<deg>(std_laguerre<deg>(std_laguerre<deg>(std_laguerre<deg>(x))))))))))));
}

template<int deg>
INLINED double aero_laguerre_12(const double x) {
	return aero_laguerre<deg>(aero_laguerre<deg>(aero_laguerre<deg>(aero_laguerre<deg>(
		aero_laguerre<deg>(aero_laguerre<deg>(aero_laguerre<deg>(aero_laguerre<deg>(
			aero_laguerre<deg>(aero_laguerre<deg>(aero_laguerre<deg>(aero_laguerre<deg>(x))))))))))));
}

template<int64_t N>
void stdh4(double* const __restrict output, const double* const __restrict input) {
	double* d_o = (double*)__builtin_assume_aligned(output, 64);
	double* d_i = (double*)__builtin_assume_aligned(input, 64);
#pragma omp parallel for
	for (int64_t i = 0; i < N; ++i) {
		d_o[i] = std_laguerre_12<12>(d_i[i]);
	}
}


template<int64_t N>
void aerobush4(double* const __restrict output, const double* const __restrict input) {
	double* d_o = (double*)__builtin_assume_aligned(output, 64);
	double* d_i = (double*)__builtin_assume_aligned(input, 64);
#pragma omp parallel for
	for (int64_t i = 0; i < N; ++i) {
		d_o[i] = aero_laguerre_12<12>(d_i[i]);
	}
}



template<int deg>
INLINED
double sin_12(const double x) {
	using AERO_SIN = aerobus::functions::sin<aerobus::i64, deg>;
	return AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(
		AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(AERO_SIN::eval(x))))))))))));
}

template<int64_t N, int deg>
INLINED
void vsin_12(const double * const __restrict in, double * const __restrict out) {
	static_assert(N % 32 == 0);
	double *d_in = (double*) __builtin_assume_aligned(in, 64);
	double *d_out = (double*) __builtin_assume_aligned(out, 64);
	#pragma omp parallel for
	for (int64_t i = 0; i < N; i += 1) {
		d_out[i] = sin_12<deg>(d_in[i]);
	}
}

INLINED
double sin_12_slow(const double x) {
	return ::sin(::sin(::sin(::sin(::sin(::sin(
		::sin(::sin(::sin(::sin(::sin(::sin(x))))))))))));
}

template<int64_t N>
INLINED
void vsin_slow(double * in, double * out) {
	static_assert(N % 32 == 0);
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
#pragma omp parallel for
	for (int i = 0; i < N; i += 1) {
		out[i] = sin_12_slow(in[i]);
	}
}


template<int64_t N>
INLINED
void fast_math_vsin_12(const double * const __restrict in, double * const __restrict out) {
	static_assert(N % 32 == 0);
	double *d_in = (double*) __builtin_assume_aligned(in, 64);
	double *d_out = (double*) __builtin_assume_aligned(out, 64);
#ifdef __AVX512F__
	#pragma omp parallel for
	for (int i = 0; i < N; i += 8)
	{
		_mm512_store_pd((d_out + i), 
			vec_sin(vec_sin(vec_sin(vec_sin(
			vec_sin(vec_sin(vec_sin(vec_sin(
			vec_sin(vec_sin(vec_sin(vec_sin(
				_mm512_load_pd((d_in + i)))))))))))))));
	}
#else
#ifdef __AVX2__
#pragma omp parallel for
	for (int i = 0; i < N; i += 4)
	{
		_mm256_store_pd((out + i),
			vec_sin(vec_sin(vec_sin(vec_sin(
			vec_sin(vec_sin(vec_sin(vec_sin(
			vec_sin(vec_sin(vec_sin(vec_sin(
				_mm256_load_pd((d_out + i)))))))))))))));
	}
#endif
#endif
}

double drand(double min, double max)
{
  double range = (max - min);
  double div = (double) RAND_MAX / range;
  return min + (rand() / div);
}

template<int deg>
void run_aero(const double* const __restrict in, double* const __restrict out_aero) {
	// warmup
	vsin_12<N, deg>(in, out_aero);
	vsin_12<N, deg>(in, out_aero);
	vsin_12<N, deg>(in, out_aero);
	double best = 1.0E9;
	for (int i = 0; i < 10; ++i) {
		auto start = std::chrono::steady_clock::now();
		vsin_12<N, deg>(in, out_aero);
		auto stop = std::chrono::steady_clock::now();
		std::chrono::duration<double> time = stop - start;
		if (time.count() < best) {
			best = time.count();
		}
	}

	printf("[aerobus deg %d] %.3e Gsin/s\n", deg, (12.0E-9 * (double) N) / best);
}

void verify(double* out_aero, double* out_std) {
	double err_aero = 0.0;
	double max_err = 0.0;
	for(size_t i = 0; i < N; ++i) {
		double err = fabs(out_aero[i] - out_std[i]);
		err_aero += err;
		if (err > max_err) {
			max_err = err;
		}
	}
	err_aero /= (double) N;
	printf("average error (vs std) : %.2e\n", err_aero);
	printf("max error (vs std) : %.2e\n", max_err);
	printf("\n");
}

int main() {
	{
		printf("#######################################################\n");
		printf("#           COMPOUND SINUS EVALUATION                 #\n");
		printf("#######################################################\n");
		printf("Evaluation of sin(x) compound 12 times with various methods\n");
		printf("\n\n");
		double* in = static_cast<double*>(aerobus::memory::aligned_malloc<double>(N, 64));
		double* out_aero = static_cast<double*>(aerobus::memory::aligned_malloc<double>(N, 64));
		double* out_std = static_cast<double*>(aerobus::memory::aligned_malloc<double>(N, 64));
		double* out_vml = static_cast<double*>(aerobus::memory::aligned_malloc<double>(N, 64));
		memset((void*)out_aero, 0, N * sizeof(double));
		memset((void*)out_std, 0, N * sizeof(double));
		memset((void*)out_vml, 0, N * sizeof(double));
#pragma omp parallel for
		for (int i = 0; i < N; ++i) {
			in[i] = drand(-0.5, 0.5); // pi / 6
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

			printf("[std math] %.3e Gsin/s\n", (12.0E-9 * (double)N) / best);
			printf("\n");
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

			printf("[std fast math] %.3e Gsin/s\n", (12.0E-9 * (double)N) / best);
			printf("\n");
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
	}
	//  laguerre
	{
		printf("#######################################################\n");
		printf("#              LAGUERRE EVALUATION                    #\n");
		printf("#######################################################\n");
		printf("Evaluation of 1-L12(x) compound 12 times with different methods\n");
		printf("Where L12 is the Laguerre polynomial of degree 12\n");
		printf("\n");
		double* input = aerobus::memory::aligned_malloc<double>(N, 1024);
		double* std_output = aerobus::memory::aligned_malloc<double>(N, 1024);
		double* aero_output = aerobus::memory::aligned_malloc<double>(N, 1024);

		for (size_t i = 0; i < N; ++i) {
			input[i] = drand(0.000000001, 0.0000001);
		}

		{
			// warmup
			stdh4<N>(std_output, input);
			stdh4<N>(std_output, input);
			stdh4<N>(std_output, input);
			std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < 4; ++i) {
				stdh4<N>(std_output, input);
			}
			std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
			::printf("time for std    : %lf\n", time_span.count());
		}

		{
			// warmup
			aerobush4<N>(aero_output, input);
			aerobush4<N>(aero_output, input);
			aerobush4<N>(aero_output, input);
			std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < 4; ++i) {
				aerobush4<N>(aero_output, input);
			}
			std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
			::printf("time for aerobus: %lf\n", time_span.count());
		}

		double max_error = 0.0;
		double average_error = 0.0;
		for (int i = 0; i < N; ++i) {
			double actual = aero_output[i];
			double expected = std_output[i];
			if (std::isnan(actual) || std::isnan(expected)) {
				printf("got NaN at %d : aero - %.15g - std - %.15g\n", i, actual, expected);
				return 1;
			}
			double error = std::abs(actual - expected);
			max_error = error > max_error ? error : max_error;
			average_error += error;
		}
		printf("Laguerre average (relative) error : %.3e\n", average_error / (double)N);
		printf("Laguerre     max (relative) error : %.3e\n", max_error);
	}
}