# Introduction

`Aerobus` is a C++-20 pure header library for general algebra on polynomials, discrete rings and associated structures.

Everything in `Aerobus` is expressed as types.

We say that again as it is the most fundamental characteristic of `Aerobus` :

***Everything is expressed as types***

The library serves two main purposes :

- Express algebra structures and associated operations in type arithmetic, compile-time;
- Provide portable and fast evaluation functions for polynomials.

It is designed to be 'quite easily' extensible.

Given these functions are "generated" at compile time and do not rely on inline assembly, they are actually platform independent, yielding exact same results if processors have same capabilities (such as Fused-Multiply-Add instructions).

## HOW TO

- Clone or download the repository somewhere, or just download `aerobus.h`
- In your code, add : `#include "aerobus.h"`
- Compile with -std=c++20 (at least) -I<install_location>
  
`Aerobus` provides a definition for low-degree (up to 997) Conway polynomials. To use them, define `AEROBUS_CONWAY_IMPORTS` before including `aerobus.h`.

### Unit Test

Install [`Cmake`](https://cmake.org/download/)
Install a recent compiler (supporting c++20), such as MSVC, G++ or Clang++

Move to the top directory then :

```bash
cmake -S . -B build
cmake --build build
cd build && ctest
```

Terminal should write :

```bash
100% tests passed, 0 tests failed out of 48
```

Alternate way :

```bash
make tests
```

From top directory.

### Benchmarks

Benchmarks are written for Intel CPUs having AVX512f and AVX512vl flags, they work only on Linux operating system using g++.

In addition of `Cmake` and compiler, install [`OpenMP`](https://www.openmp.org/resources/openmp-compilers-tools/).
Then move to top directory :

```bash
rm -rf build 
mkdir build
cd build
cmake ..
make aerobus_benchmarks
./aerobus_benchmarks
```

## Structures

### Predefined discrete euclidean domains

`Aerobus` predefines several simple euclidean domains, such as :

- `aerobus::i32` : integers (32 bits)
- `aerobus::i64` : integers (64 bits)
- `aerobus::zpz<p>` : integers modulo p (prime number) on 32 bits

All these types represent the Ring, meaning the algebraic structure. They have a nested type `val<i>` where `i` is a scalar native value (int32_t or int64_t) to represent actual values in the ring.
They have the following "operations", required by the IsEuclideanDomain concept :

- `add_t` : a type (specialization of val), representing addition between two values
- `sub_t` : a type (specialization of val), representing subtraction between two values
- `mul_t` : a type (specialization of val), representing multiplication between two values
- `div_t` : a type (specialization of val), representing division between two values
- `mod_t` : a type (specialization of val), representing modulus between two values

and the following "elements" :

- one : the neutral element for multiplication, val<1>
- zero : the neutral element for addition, val<0>

### Polynomials

`Aerobus` defines polynomials as a variadic template structure, with coefficient in an arbitrary discrete euclidean domain. As `i32` or `i64`, they are given same operations and elements, which make them a euclidean domain by themselves. Similarly, `aerobus::polynomial` represents the algebraic structure, actual values are in `aerobus::polynomial::val`.

In addition, values have an evaluation function :

```cpp
template<typename valueRing> static constexpr valueRing eval(const valueRing& x) {...}
```

Which can be used at compile time (constexpr evaluation) or runtime.

### Known polynomials

`Aerobus` predefines some well known families of polynomials, such as Hermite or Bernstein :

```cpp
using B23 = aerobus::known_polynomials::bernstein<2, 3>; // 3X^2(1-X)
constexpr float x = B32::eval(2.0F); // -12
```

They have their coefficients either in `aerobus::i64` or `aerobus::q64`.
Complete list is (but is meant to be extended):

- `chebyshev_T`
- `chebyshev_U`
- `laguerre`
- `hermite_prob`
- `hermite_phys`
- `bernstein`
- `legendre`
- `bernoulli`

### Conway polynomials

When the tag `AEROBUS_CONWAY_IMPORTS` is defined at compile time (`-DAEROBUS_CONWAY_IMPORTS`), `aerobus` provides definition for all Conway polynomials `CP(p, n)` for `p` up to 997 and low values for `n` (usually less than 10).

They can be used to construct finite fields of order $p^n$ ($\mathbb{F}_{p^n}$):

```cpp
using F2 = zpz<2>;
using PF2 = polynomial<F2>;
using F4 = Quotient<PF2, ConwayPolynomial<2, 2>::type>;
```

### Taylor series

`Aerobus` provides definition for Taylor expansion of known functions. They are all templates in two parameters, degree of expansion (`size_t`) and Integers (`typename`). Coefficients then live in `FractionField<Integers>`.

They can be used and evaluated:

```cpp
using namespace aerobus;
using aero_atanh = atanh<i64, 6>;
constexpr float val = aero_atanh::eval(0.1F); // approximation of arctanh(0.1) using taylor expansion of degree 6
```

Exposed functions are:

- `exp`
- `expm1` $e^x - 1$
- `lnp1` $\ln(x+1)$
- `geom` $\frac{1}{1-x}$
- `sin`
- `cos`
- `tan`
- `sh`
- `cosh`
- `tanh`
- `asin`
- `acos`
- `acosh`
- `asinh`
- `atanh`

Having the capacity of specifying the degree is very important, as users may use other formats than `float64` or `float32` which require higher or lower degree to achieve correct or acceptable precision.

It's possible to define Taylor expansion by implementing a `coeff_at` structure which must meet the following requirement :

- Being template in Integers (`typename`) and index (`size_t`);
- Exposing a type alias `type`, some specialization of `FractionField<Integers>::val`.

For example, to define the serie $1+x+x^2+x^3+\ldots$, users may write:

```cpp
template<typename Integers, size_t i>
struct my_coeff_at {
    using type = typename FractionField<Integers>::one;
};

template<typename Integers, size_t degree>
using my_serie = taylor<Integers, my_coeff_at, degree>;

static constexpr double x = my_serie<i64, 3>::eval(3.0);
```

On x86-64 and CUDA platforms at least, using proper compiler directives, these functions yield very performant assembly, similar or better than standard library implementation in fast math. For example, this code:

```cpp
double compute_expm1(const size_t N, double* in, double* out) {
    using V = aerobus::expm1<aerobus::i64, 13>;
    for (size_t i = 0; i < N; ++i) {
        out[i] = V::eval(in[i]);
    }
}
```

Yields this assembly (clang 17, `-mavx2 -O3`) where we can see a pile of Fused-Multiply-Add vector instructions, generated because we unrolled completely the Horner evaluation loop:

```nasm
compute_expm1(unsigned long, double const*, double*):
  lea     rax, [rdi-1]
  cmp     rax, 2
  jbe     .L5
  mov     rcx, rdi
  xor     eax, eax
  vxorpd  xmm1, xmm1, xmm1
  vbroadcastsd    ymm14, QWORD PTR .LC1[rip]
  vbroadcastsd    ymm13, QWORD PTR .LC3[rip]
  shr     rcx, 2
  vbroadcastsd    ymm12, QWORD PTR .LC5[rip]
  vbroadcastsd    ymm11, QWORD PTR .LC7[rip]
  sal     rcx, 5
  vbroadcastsd    ymm10, QWORD PTR .LC9[rip]
  vbroadcastsd    ymm9, QWORD PTR .LC11[rip]
  vbroadcastsd    ymm8, QWORD PTR .LC13[rip]
  vbroadcastsd    ymm7, QWORD PTR .LC15[rip]
  vbroadcastsd    ymm6, QWORD PTR .LC17[rip]
  vbroadcastsd    ymm5, QWORD PTR .LC19[rip]
  vbroadcastsd    ymm4, QWORD PTR .LC21[rip]
  vbroadcastsd    ymm3, QWORD PTR .LC23[rip]
  vbroadcastsd    ymm2, QWORD PTR .LC25[rip]
.L3:
  vmovupd ymm15, YMMWORD PTR [rsi+rax]
  vmovapd ymm0, ymm15
  vfmadd132pd     ymm0, ymm14, ymm1
  vfmadd132pd     ymm0, ymm13, ymm15
  vfmadd132pd     ymm0, ymm12, ymm15
  vfmadd132pd     ymm0, ymm11, ymm15
  vfmadd132pd     ymm0, ymm10, ymm15
  vfmadd132pd     ymm0, ymm9, ymm15
  vfmadd132pd     ymm0, ymm8, ymm15
  vfmadd132pd     ymm0, ymm7, ymm15
  vfmadd132pd     ymm0, ymm6, ymm15
  vfmadd132pd     ymm0, ymm5, ymm15
  vfmadd132pd     ymm0, ymm4, ymm15
  vfmadd132pd     ymm0, ymm3, ymm15
  vfmadd132pd     ymm0, ymm2, ymm15
  vfmadd132pd     ymm0, ymm1, ymm15
  vmovupd YMMWORD PTR [rdx+rax], ymm0
  add     rax, 32
  cmp     rcx, rax
  jne     .L3
  mov     rax, rdi
  and     rax, -4
  vzeroupper
```

## Operations

### Field of fractions

Given a set (type) satisfies the `IsEuclideanDomain` concept, `Aerobus` allows to define its [field of fractions](https://en.wikipedia.org/wiki/Field_of_fractions).

This new type is again a euclidean domain, especially a field, and therefore we can define polynomials over it.

For example, integers modulo `p` is not a field when `p` is not prime. We then can define its field of fraction and polynomials over it this way:

```cpp
using namespace aerobus;
using ZmZ = zpz<8>;
using Fzmz = FractionField<ZmZ>;
using Pfzmz = polynomial<Fzmz>;
```

The same operation would stand for any set that users would have implemented in place of `ZmZ`.

For example, we can easily define [rational functions](https://en.wikipedia.org/wiki/Rational_function) by taking the ring of fractions of polynomials:

```cpp
using namespace aerobus;
using RF64 = FractionField<polynomial<q64>>;
```

Which also have an evaluation function, as polynomial do.

### Quotient

Given a ring `R`, `Aerobus` provides automatic implementation for [quotient ring](https://en.wikipedia.org/wiki/Quotient_ring) $R/X$ where X is a principal ideal generated by some element, as we know this kind of ideal is two-sided as long as `R` is commutative (and we assume it is).

For example, if we want `R` to be $\mathbb{Z}$ represented as `aerobus::i64`, we can express arithmetic modulo 17 using:

```cpp
using namespace aerobus;
using ZpZ = Quotient<i64, i64::val<17>>;
```

As we could have using `zpz<17>`.

This is mainly used to define finite fields of order $p^n$ using Conway polynomials but may have other applications.

## Misc

### Continued Fractions

`Aerobus` gives an implementation for [continued fractions](https://en.wikipedia.org/wiki/Continued_fraction). It can be used this way:

```cpp
using namespace aerobus;
using T = ContinuedFraction<1,2,3,4>;
constexpr double x = T::val;
```

As practical examples, `aerobus` gives continued fractions of $\pi$, $e$, $\sqrt{2}$ and $\sqrt{3}$:

```cpp
constexpr double A_SQRT3 = aerobus::SQRT3_fraction::val; // 1.7320508075688772935
```

## CUDA

When compiled with `nvcc` and the flag `WITH_CUDA_FP16`, `Aerobus` provides some kind of support of 16 bits integers and floats (aka `__half`).

Unfortunately, NVIDIA did not put enough constexpr in its `cuda_fp16.h` header, so we had to implement our own constexpr static_cast from int16_t to `__half` to make integers polynomials work with `__half`. See [this bug](https://developer.nvidia.com/bugs/4863696).

More, it's (at this time), not possible to make it work for `__half2` because of [another bug](https://developer.nvidia.com/bugs/4872028).

Please push to make these bug fixed by NVIDIA.

[![DOI](https://zenodo.org/badge/499577459.svg)](https://zenodo.org/badge/latestdoi/499577459)
