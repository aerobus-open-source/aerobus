<<<<<<< HEAD
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include "lib.h"

INLINED
double expm1_12(const double x) {
	using V = aerobus::expm1<aerobus::i64, 13>;
	return V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(
		V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(x))))))))))));
}

INLINED
void vexpm1_12(double * in, double * out, size_t N) {
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		out[i] = expm1_12(in[i]);
	}
}

int main() {
	const size_t N = 100000000;
	double* in = static_cast<double*>(aligned_alloc(64, N * sizeof(double)));
	double* out = static_cast<double*>(aligned_alloc(64, N * sizeof(double)));
	for(int i = 0; i < N; ++i) {
		in[i] = 0.1;
	}
	// warmup
	vexpm1_12(in, out, N);
	vexpm1_12(in, out, N);
	vexpm1_12(in, out, N);
	double best = 1.0E9;
	for (int i = 0; i < 10; ++i) {
		auto start = std::chrono::steady_clock::now();
		vexpm1_12(in, out, N);
		auto stop = std::chrono::steady_clock::now();
		std::chrono::duration<double> time = stop - start;
		if (time.count() < best) {
			best = time.count();
		}
	}

	printf("time for 100E6 exp12 : %lf\n", best);
	printf("achieved Gflops : %lf\n", ((double)N * 13.0 * 12.0E-9 / best));
=======
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include "lib.h"


INLINED
double expm1_12(const double x) {
	using V = aerobus::expm1<aerobus::i64, 13>;
	return V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(
		V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(x))))))))))));
}

INLINED
void vexpm1_12(double * in, double * out, size_t N) {
	in = (double*) __builtin_assume_aligned(in, 64);
	out = (double*) __builtin_assume_aligned(out, 64);
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		out[i] = expm1_12(in[i]);
	}
}

int main() {
	const size_t N = 100000000;
	double* in = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	double* out = static_cast<double*>(aerobus::aligned_malloc<double>(N, 64));
	for(int i = 0; i < N; ++i) {
		in[i] = 0.1;
	}
	// warmup
	vexpm1_12(in, out, N);
	vexpm1_12(in, out, N);
	vexpm1_12(in, out, N);
	double best = 1.0E9;
	for (int i = 0; i < 10; ++i) {
		auto start = std::chrono::steady_clock::now();
		vexpm1_12(in, out, N);
		auto stop = std::chrono::steady_clock::now();
		std::chrono::duration<double> time = stop - start;
		if (time.count() < best) {
			best = time.count();
		}
	}

	printf("time for 100E6 exp12 : %lf\n", best);
	printf("achieved Gflops : %lf\n", ((double)N * 13.0 * 12.0E-9 / best));
>>>>>>> d26c25a (merge headers back into one big file)
}