# aerobus

Polynomials, rings and finite fields computed by the C++ compiler.

`aerobus` is a header-only C++20 library where every algebraic value is a type: an integer modulo p, a rational number, a polynomial, an element of a finite field. The compiler does the arithmetic, so results are exact and can be checked with `static_assert`. What reaches the program is constants: lookup tables, or polynomials evaluated by a fully unrolled Horner scheme on `double`, `float` or CUDA `__half`.

Version 2.0 breaks the API of 1.x: see [CHANGELOG.md](CHANGELOG.md).

## Three examples

### Finite fields, for error-correcting codes

```cpp
using F2 = aerobus::zpz<2>;
using GF256 = aerobus::GaloisField<2, 8>;  // F2[x] / (x^8 + x^4 + x^3 + x^2 + 1)
using x = GF256::inject_ring_t<aerobus::polynomial<F2>::X>;

// x^8 = x^4 + x^3 + x^2 + 1: the byte 0x1D of the QR code standard
static_assert(std::is_same_v<aerobus::pow_t<x, 8>,
                             GF256::inject_ring_t<aerobus::make_int_polynomial_t<F2, 1, 1, 1, 0, 1>>>);
static_assert(std::is_same_v<aerobus::pow_t<x, 255>, GF256::one>);
```

`GaloisField<p, n>` is the field with pⁿ elements, built from the Conway polynomial CP(p, n). For p = 2 and n = 8 it is the field of QR codes and RAID-6. The article [Reed-Solomon with a field the compiler builds and proves](articles/reed_solomon_qr.md) builds a QR code error-correction codec on it, with tables computed and proven at compile time, and compares it with zxing-cpp and Schifra. [Run this example on Compiler Explorer](https://godbolt.org/clientstate/eyJzZXNzaW9ucyI6W3siaWQiOjEsImxhbmd1YWdlIjoiYysrIiwic291cmNlIjoiLy8gR0YoMjU2KSwgdGhlIGZpZWxkIG9mIFFSIGNvZGVzIGFuZCBSQUlELTYsIGNvbXB1dGVkIGJ5IHRoZSBjb21waWxlciB3aXRoIGFlcm9idXNcbiNpbmNsdWRlIDxjc3RkaW8%2BXG4jZGVmaW5lIEFFUk9CVVNfQ09OV0FZX0lNUE9SVFNcbiNpbmNsdWRlIDxodHRwczovL3Jhdy5naXRodWJ1c2VyY29udGVudC5jb20vYWVyb2J1cy1vcGVuLXNvdXJjZS9hZXJvYnVzL3YyLjAuMC9zcmMvYWVyb2J1cy5oPlxuXG51c2luZyBGMiA9IGFlcm9idXM6OnpwejwyPjtcbnVzaW5nIEdGMjU2ID0gYWVyb2J1czo6R2Fsb2lzRmllbGQ8MiwgOD47ICAvLyBGMlt4XSAvICh4XjggKyB4XjQgKyB4XjMgKyB4XjIgKyAxKVxudXNpbmcgeCA9IEdGMjU2OjppbmplY3RfcmluZ190PGFlcm9idXM6OnBvbHlub21pYWw8RjI%2BOjpYPjtcblxuLy8geF44ID0geF40ICsgeF4zICsgeF4yICsgMSwgdGhlIGJ5dGUgMHgxRCBvZiB0aGUgUVIgY29kZSBzdGFuZGFyZFxuc3RhdGljX2Fzc2VydChzdGQ6OmlzX3NhbWVfdjxhZXJvYnVzOjpwb3dfdDx4LCA4PixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgR0YyNTY6OmluamVjdF9yaW5nX3Q8YWVyb2J1czo6bWFrZV9pbnRfcG9seW5vbWlhbF90PEYyLCAxLCAxLCAxLCAwLCAxPj4%2BKTtcbi8vIHggaGFzIG9yZGVyIDI1NTogaXQgZ2VuZXJhdGVzIHRoZSAyNTUgbm9uemVybyBlbGVtZW50c1xuc3RhdGljX2Fzc2VydChzdGQ6OmlzX3NhbWVfdjxhZXJvYnVzOjpwb3dfdDx4LCAyNTU%2BLCBHRjI1Njo6b25lPik7XG5zdGF0aWNfYXNzZXJ0KCFzdGQ6OmlzX3NhbWVfdjxhZXJvYnVzOjpwb3dfdDx4LCA4NT4sIEdGMjU2OjpvbmU%2BICYmICFzdGQ6OmlzX3NhbWVfdjxhZXJvYnVzOjpwb3dfdDx4LCA1MT4sIEdGMjU2OjpvbmU%2BXG4gICAgICAgICAgICAgICYmICFzdGQ6OmlzX3NhbWVfdjxhZXJvYnVzOjpwb3dfdDx4LCAxNT4sIEdGMjU2OjpvbmU%2BKTtcblxuaW50IG1haW4oKSB7XG4gICAgc3RkOjpwcmludGYoXCJ4XjggPSAlcyBpbiBHRigyNTYpXFxuXCIsIGFlcm9idXM6OnBvd190PHgsIDg%2BOjp0eXBlOjp0b19zdHJpbmcoKS5jX3N0cigpKTtcbn1cbiIsImNvbXBpbGVycyI6W3siaWQiOiJnMTUyIiwib3B0aW9ucyI6Ii1zdGQ9YysrMjAgLU8yIC1mdGVtcGxhdGUtZGVwdGg9MzAwMDAifV0sImV4ZWN1dG9ycyI6W3siY29tcGlsZXIiOnsiaWQiOiJnMTUyIiwib3B0aW9ucyI6Ii1zdGQ9YysrMjAgLU8yIC1mdGVtcGxhdGUtZGVwdGg9MzAwMDAiLCJsaWJzIjpbXX0sImFyZ3VtZW50cyI6IiIsInN0ZGluIjoiIn1dfV19).

### Orthogonal polynomials with exact coefficients

```cpp
using P6 = aerobus::known_polynomials::legendre<6>;  // (231x^6 - 315x^4 + 105x^2 - 5) / 16
static_assert(std::is_same_v<P6::coeff_at_t<6>, aerobus::make_q64_t<231, 16>>);

double p6(double x) { return P6::eval(x); }
```

Coefficients come from the recurrence in exact rational arithmetic and are rounded once, when converted to `double`. `p6` compiles to constants and multiply-adds, with no loop. Chebyshev, Hermite, Laguerre, Bernstein, Bernoulli and others are predefined. [Run this example on Compiler Explorer](https://godbolt.org/clientstate/eyJzZXNzaW9ucyI6W3siaWQiOjEsImxhbmd1YWdlIjoiYysrIiwic291cmNlIjoiLy8gTGVnZW5kcmUgcG9seW5vbWlhbCBQNiB3aXRoIGV4YWN0IHJhdGlvbmFsIGNvZWZmaWNpZW50cywgZXZhbHVhdGVkIGFzIGFuIHVucm9sbGVkIEhvcm5lciBjaGFpblxuI2luY2x1ZGUgPGNtYXRoPlxuI2luY2x1ZGUgPGNzdGRpbz5cbiNpbmNsdWRlIDxodHRwczovL3Jhdy5naXRodWJ1c2VyY29udGVudC5jb20vYWVyb2J1cy1vcGVuLXNvdXJjZS9hZXJvYnVzL3YyLjAuMC9zcmMvYWVyb2J1cy5oPlxuXG51c2luZyBQNiA9IGFlcm9idXM6Omtub3duX3BvbHlub21pYWxzOjpsZWdlbmRyZTw2PjsgIC8vICgyMzF4XjYgLSAzMTV4XjQgKyAxMDV4XjIgLSA1KSAvIDE2XG5zdGF0aWNfYXNzZXJ0KHN0ZDo6aXNfc2FtZV92PFA2Ojpjb2VmZl9hdF90PDY%2BLCBhZXJvYnVzOjptYWtlX3E2NF90PDIzMSwgMTY%2BPik7XG5zdGF0aWNfYXNzZXJ0KHN0ZDo6aXNfc2FtZV92PFA2Ojpjb2VmZl9hdF90PDA%2BLCBhZXJvYnVzOjptYWtlX3E2NF90PC01LCAxNj4%2BKTtcblxuZG91YmxlIHA2KGRvdWJsZSB4KSB7IHJldHVybiBQNjo6ZXZhbCh4KTsgfSAgLy8gc2VlIHRoZSBhc3NlbWJseTogY29uc3RhbnRzIGFuZCBtdWx0aXBseS1hZGRzIG9ubHlcblxuaW50IG1haW4oKSB7XG4gICAgc3RkOjpwcmludGYoXCJQNigwLjMpID0gJS4xN2dcXG5zdGQ6OmxlZ2VuZHJlKDYsIDAuMykgPSAlLjE3Z1xcblwiLCBwNigwLjMpLCBzdGQ6OmxlZ2VuZHJlKDYsIDAuMykpO1xufVxuIiwiY29tcGlsZXJzIjpbeyJpZCI6ImcxNTIiLCJvcHRpb25zIjoiLXN0ZD1jKysyMCAtTzIgLWZ0ZW1wbGF0ZS1kZXB0aD0zMDAwMCJ9XSwiZXhlY3V0b3JzIjpbeyJjb21waWxlciI6eyJpZCI6ImcxNTIiLCJvcHRpb25zIjoiLXN0ZD1jKysyMCAtTzIgLWZ0ZW1wbGF0ZS1kZXB0aD0zMDAwMCIsImxpYnMiOltdfSwiYXJndW1lbnRzIjoiIiwic3RkaW4iOiIifV19XX0%3D).

### Compensated evaluation in `float`

```cpp
// (x - 1)^6, expanded: catastrophic cancellation near x = 1
using P = aerobus::make_int_polynomial_t<aerobus::i64, 1, -6, 15, -20, 15, -6, 1>;

float naive(float x) { return P::eval(x); }
float accurate(float x) { return P::compensated_eval(x); }
```

At x = 1.05F the exact value is 1.562491e-8. Horner in `float` returns -5.96e-7 with GCC and 8.56e-7 with clang; `compensated_eval` returns 1.562495e-8. The compensated Horner scheme is about as accurate as Horner in twice the working precision, which matters on GPUs whose `double` throughput is a small fraction of `float`. [Run this example on Compiler Explorer](https://godbolt.org/clientstate/eyJzZXNzaW9ucyI6W3siaWQiOjEsImxhbmd1YWdlIjoiYysrIiwic291cmNlIjoiLy8gQ29tcGVuc2F0ZWQgSG9ybmVyOiBmbG9hdCBldmFsdWF0aW9uIGFib3V0IGFzIGFjY3VyYXRlIGFzIEhvcm5lciBpbiBkb3VibGVcbiNpbmNsdWRlIDxjc3RkaW8%2BXG4jaW5jbHVkZSA8aHR0cHM6Ly9yYXcuZ2l0aHVidXNlcmNvbnRlbnQuY29tL2Flcm9idXMtb3Blbi1zb3VyY2UvYWVyb2J1cy92Mi4wLjAvc3JjL2Flcm9idXMuaD5cblxuLy8gKHggLSAxKV42LCBleHBhbmRlZDogY2F0YXN0cm9waGljIGNhbmNlbGxhdGlvbiBuZWFyIHggPSAxXG51c2luZyBQID0gYWVyb2J1czo6bWFrZV9pbnRfcG9seW5vbWlhbF90PGFlcm9idXM6Omk2NCwgMSwgLTYsIDE1LCAtMjAsIDE1LCAtNiwgMT47XG5cbmZsb2F0IG5haXZlKGZsb2F0IHgpIHsgcmV0dXJuIFA6OmV2YWwoeCk7IH1cbmZsb2F0IGFjY3VyYXRlKGZsb2F0IHgpIHsgcmV0dXJuIFA6OmNvbXBlbnNhdGVkX2V2YWwoeCk7IH1cblxuaW50IG1haW4oKSB7XG4gICAgY29uc3QgZmxvYXQgeCA9IDEuMDVGO1xuICAgIGNvbnN0IGRvdWJsZSBkID0gc3RhdGljX2Nhc3Q8ZG91YmxlPih4KSAtIDE7ICAvLyBleGFjdDogeCAtIDEgaXMgcmVwcmVzZW50YWJsZVxuICAgIHN0ZDo6cHJpbnRmKFwiZXhhY3QgICAgICAgICUuNmVcXG5Ib3JuZXIgICAgICAgJS42ZVxcbmNvbXBlbnNhdGVkICAlLjZlXFxuXCIsIGQgKiBkICogZCAqIGQgKiBkICogZCxcbiAgICAgICAgICAgICAgICBzdGF0aWNfY2FzdDxkb3VibGU%2BKG5haXZlKHgpKSwgc3RhdGljX2Nhc3Q8ZG91YmxlPihhY2N1cmF0ZSh4KSkpO1xufVxuIiwiY29tcGlsZXJzIjpbeyJpZCI6ImcxNTIiLCJvcHRpb25zIjoiLXN0ZD1jKysyMCAtTzIgLWZ0ZW1wbGF0ZS1kZXB0aD0zMDAwMCJ9XSwiZXhlY3V0b3JzIjpbeyJjb21waWxlciI6eyJpZCI6ImcxNTIiLCJvcHRpb25zIjoiLXN0ZD1jKysyMCAtTzIgLWZ0ZW1wbGF0ZS1kZXB0aD0zMDAwMCIsImxpYnMiOltdfSwiYXJndW1lbnRzIjoiIiwic3RkaW4iOiIifV19XX0%3D).

## Getting started

- Clone or download the repository somewhere, or just download `src/aerobus.h`
- In your code, add: `#include "aerobus.h"`
- Compile with `-std=c++20` (at least) and `-I<install_location>`
- Deep compile-time computations (high degree Taylor series, large Conway polynomials) may need a larger template depth than the compiler default: `-ftemplate-depth=30000` with GCC and clang.

With CMake, either add the repository with `add_subdirectory`, or install it and use the package:

```cmake
find_package(Aerobus 2 REQUIRED)
target_link_libraries(my_target PRIVATE aerobus::aerobus)
```

`Aerobus` provides the Conway polynomials `CP(p, n)` for primes `p` below 1000 (degrees up to 20 for small `p`, 9 for the largest). To use them, define `AEROBUS_CONWAY_IMPORTS` before including `aerobus.h`. The table is generated by `tools/generate_conway.py` from `imports/conwaypolynomials.txt`.

Supported compilers: GCC 13 or newer, clang 17 or newer, and nvcc 12.6 or newer for the CUDA code paths (tested with GCC 15.2, clang 21.1 and nvcc 13.2).

### Unit tests

Install [`CMake`](https://cmake.org/download/) (3.14 or newer) and a C++20 compiler, then from the top directory:

```bash
cmake -S . -B build
cmake --build build
cd build && ctest
```

All tests should pass. The tests include the examples (each checks its own result), compile-fail tests for the static assertions, a check that every C++ block of this README compiles, and, when nvcc is found, host-side checks of the CUDA fp16 code.

Alternate way, which runs the tests with clang then GCC:

```bash
make tests
```

### Benchmarks

The benchmarks use Google's [Benchmark library](https://github.com/google/benchmark) (found on the system or downloaded) and [OpenMP](https://www.openmp.org/resources/openmp-compilers-tools/) when available. They are compiled with `-march=native` (set `AEROBUS_BENCHMARK_FLAGS` to change it):

```bash
cmake -S . -B build -DAEROBUS_BUILD_BENCHMARKS=ON
cmake --build build --target benchmarks
./build/benchmarks
```

## Limitations

- **64-bit coefficients.** Coefficients are computed exactly in `i64` (or `i32`, `i16`), so degrees are bounded, unless the integers are `bigint`, at a compile time cost of about 2 to 4 times. Each predefined family requires `aerobus::coefficients_fit`: a degree above the measured limit is a constraint failure that names it, and `requires { typename legendre<22>; }` is false. The limits in `i64` are 52 for `chebyshev_T`, 51 for `chebyshev_U`, 43 for `bernstein`, 34 for `bernoulli`, 32 for `hermite_prob`, 25 for `hermite_phys` and `touchard`, 21 for `legendre`, 17 for `bessel`, 14 for `laguerre`, 21 for the Taylor series of `cos` and `cosh` and 20 for the other Taylor series; `lnp1`, `atan`, `atanh` and `geometric_sum` have no limit. The limits in `i32` and `i16` are in the `max_degree` table of `src/aerobus.h`.
- **Taylor series are not a libm.** They have no argument reduction: accurate near 0 only. Use them as building blocks, not as a replacement for `std::exp`.
- **Results depend on FMA contraction.** Whether `x * y + z` becomes a fused multiply-add depends on the target and the compiler flags, so the last bits can differ between compilers and options (see the example above).
- **Compile time.** Deep computations cost seconds and may need `-ftemplate-depth=30000`. The Reed-Solomon example compiles in about 4 s with GCC and 7 s with clang; MSVC builds the library but does not finish that example.
- **Not constant-time.** Table lookups depend on the data and leak through the cache: fine for error correction, not for cryptography.
- **`compensated_eval` needs a correctly rounded `std::fma`**, which is slow on processors without a fused multiply-add instruction.

The rest of this page is the reference, structure by structure.

## Structures

### Predefined discrete euclidean domains

`Aerobus` predefines several simple euclidean domains, such as:

- `aerobus::i16`, `aerobus::i32`, `aerobus::i64`: integers on 16, 32 and 64 bits (`aerobus::integers<Int>`)
- `aerobus::zpz<p>`: integers modulo `p`, on 32 bits; a field when `p` is prime

All these types represent the Ring, meaning the algebraic structure. They have a nested type `val<i>` where `i` is a scalar native value to represent actual values in the ring.
They have the following "operations", required by the `IsEuclideanDomain` concept:

- `add_t`: a type (specialization of val), representing addition between two values
- `sub_t`: a type (specialization of val), representing subtraction between two values
- `mul_t`: a type (specialization of val), representing multiplication between two values
- `div_t`: a type (specialization of val), representing division between two values (by the inverse in `zpz<p>`)
- `mod_t`: a type (specialization of val), representing modulus between two values

and the following "elements":

- `one`: the neutral element for multiplication, `val<1>`
- `zero`: the neutral element for addition, `val<0>`

Every value has a single type: `zpz<7>::val<-1>` is `zpz<7>::val<6>`. The value-level helpers infer the ring from their arguments (`aerobus::add_t<a, b>`, `aerobus::divide_t<a, b>`, `aerobus::gcd_t<a, b>`, `aerobus::pow_t<a, n>`):

```cpp
using namespace aerobus;
using Z7 = zpz<7>;
static_assert(std::is_same_v<Z7::val<-1>, Z7::val<6>>);
static_assert(std::is_same_v<divide_t<Z7::one, Z7::val<2>>, Z7::val<4>>);  // 2 * 4 = 8 = 1
static_assert(gcd_t<i32::val<12>, i32::val<18>>::v == 6);
static_assert(pow_t<i64::val<2>, 10>::v == 1024);
```

### Arbitrary precision integers

`aerobus::bigint` is an euclidean domain of integers without bound, usable wherever `i64` is: as polynomial coefficients, in `FractionField<bigint>`, as the integers of the predefined families and Taylor series, which are then not limited in degree:

```cpp
using aerobus::bigint;
using N = bigint::from_string_t<"-123456789012345678901234567890">;
using M = bigint::mul_t<N, bigint::from_hex_t<"FFFFFFFFFFFFFFFF">>;
static_assert(bigint::eq_v<bigint::div_t<M, N>, bigint::from_hex_t<"FFFFFFFFFFFFFFFF">>);
static_assert(bigint::mod_t<bigint::inject_constant_t<-7>, bigint::inject_constant_t<2>>::get<int>() == -1);
using P30 = aerobus::known_polynomials::legendre<30, bigint>;  // exact coefficients, beyond the i64 limit of 21
```

Values are `bigint::val<negative, digits...>`, digits in base 2^32, most significant first. The arithmetic runs in `constexpr` functions on arrays of digits rather than through recursive templates, and division is Knuth's algorithm D, so dividing a 40-digit number by a 20-digit one adds about 0.04 s to the compilation. As for `i64`, division truncates toward zero and the remainder has the sign of the dividend. `to_string()` and `to_hex()` print values, `get<T>()` converts them.

### Polynomials

`Aerobus` defines polynomials as a variadic template structure, with coefficients in an arbitrary discrete euclidean domain. As `i32` or `i64`, they are given the same operations and elements, which make them a euclidean domain by themselves. Similarly, `aerobus::polynomial<Ring>` represents the algebraic structure, actual values are in `aerobus::polynomial<Ring>::val`, which takes the coefficients highest degree first and removes leading zeros.

Division works over fields; over other rings (such as `i64`) it works when every step is exact, and stops with a `static_assert` otherwise.

Values have an evaluation function, usable at compile time (constexpr evaluation) or at run time:

```cpp
using P = aerobus::make_int_polynomial_t<aerobus::i64, 3, 2, 1>;  // 3x^2 + 2x + 1
constexpr double y = P::eval(2.0);  // 17
static_assert(y == 17.0);
```

`eval<T>` converts the coefficients to `T`: for coefficients in `zpz<p>` this uses their representatives in `[0, p)`, which is not arithmetic in `Z/pZ`. To evaluate inside the ring, use `value_at_t`:

```cpp
using F2 = aerobus::zpz<2>;
using P = aerobus::make_int_polynomial_t<F2, 1, 0, 0, 1>;  // x^3 + 1
static_assert(P::value_at_t<F2::one>::v == 0);  // 1 + 1 = 0 in Z/2Z
```

`compensated_eval` evaluates with the compensated Horner scheme (as accurate as Horner in twice the working precision), at run time and for floating point types. It relies on a correctly rounded fused multiply-add (`std::fma`), which is slow on processors without one.

### Known polynomials

`Aerobus` predefines some well known families of polynomials, such as Hermite or Bernstein:

```cpp
using B23 = aerobus::known_polynomials::bernstein<2, 3>;  // 3X^2(1-X)
constexpr float x = B23::eval(2.0F);  // -12
static_assert(x == -12.0F);
```

They have their coefficients either in `aerobus::i64` or `aerobus::q64` (the integer ring is a template parameter).
Complete list is (but is meant to be extended):

- `chebyshev_T`, `chebyshev_U`
- `laguerre`
- `hermite_prob`, `hermite_phys`
- `bernstein`
- `legendre`
- `bernoulli`
- `allone`
- `bessel`
- `touchard`
- `abel`

### Conway polynomials

When `AEROBUS_CONWAY_IMPORTS` is defined at compile time (`-DAEROBUS_CONWAY_IMPORTS`), `aerobus` provides all Conway polynomials `CP(p, n)` for `p` below 1000 and low values of `n`.

They can be used to construct finite fields of order $p^n$ ($\mathbb{F}_{p^n}$):

```cpp
using F2 = aerobus::zpz<2>;
using PF2 = aerobus::polynomial<F2>;
using F4 = aerobus::Quotient<PF2, aerobus::ConwayPolynomial<2, 2>::type, true>;  // true: X is irreducible
using F9 = aerobus::GaloisField<3, 2>;  // the same construction, for p = 3 and n = 2
static_assert(F4::is_field && F9::is_field);
```

### Factorization over finite fields

Polynomials with coefficients in a finite field can be tested for irreducibility and factored, at compile time. The finite fields are `zpz<p>` with `p` prime, the quotients of integers by a prime, and `Quotient<polynomial<zpz<p>>, X, true>`, which includes `GaloisField<p, n>` (concept `IsFiniteField`):

```cpp
using F2 = aerobus::zpz<2>;
using P = aerobus::make_int_polynomial_t<F2, 1, 0, 0, 0, 0, 0, 0, 1>;  // x^7 + 1
using F = aerobus::factorization_t<P>;  // (x + 1)(x^3 + x + 1)(x^3 + x^2 + 1)
static_assert(F::factors::length == 3);
static_assert(std::is_same_v<F::type, P>);  // the product of the factors
static_assert(aerobus::is_irreducible_v<aerobus::ConwayPolynomial<2, 8>::type>);
```

- `is_irreducible_v<P>`: Rabin's test.
- `squarefree_t<P>`: squarefree decomposition over any field, including `q64` (Musser's algorithm, with p-th roots in characteristic p).
- `factorization_t<P>`: squarefree decomposition, distinct degree factorization, then Cantor-Zassenhaus with a deterministic sequence of candidates instead of random ones.

The result is a `factorization<unit, type_list<factor<F, e>...>>`. `unit` is the leading coefficient, the factors are monic and sorted, so that equal factorizations are the same type, and the member `type` is their product. `field_order_v`, `characteristic_v` and `field_element_t<F, i>` describe and enumerate the field. Exponents such as x^(q^d) are reached by repeated q-th powers, never computed as integers. The cost is compile time: factoring x^15 - 1 over GF(4) takes a few seconds and about 1 GB of compiler memory.

The example [`examples/bch_qr_format.cpp`](examples/bch_qr_format.cpp) factors x^15 - 1 over F2, multiplies the minimal polynomials of alpha .. alpha^6 into the generator of the BCH(15, 5) code that protects the format information of QR codes (0x537), and rebuilds the 32 format words of the standard. The article [The QR format code is a factor of x¹⁵ − 1](articles/bch_qr_format.md) walks through it.

### Taylor series

`Aerobus` provides Taylor expansions of known functions. They are all templates in two parameters: the integers (`typename`) and the degree of expansion (`size_t`). Coefficients live in `FractionField<Integers>`.

They can be used and evaluated:

```cpp
using aero_atanh = aerobus::atanh<aerobus::i64, 6>;
constexpr float val = aero_atanh::eval(0.1F);  // approximation of arctanh(0.1) using taylor expansion of degree 6
```

Qualify them with `aerobus::`: after `using namespace aerobus`, names such as `exp`, `sin` or `atanh` are ambiguous with the C functions of `<math.h>`.

Exposed functions are:

- `exp`
- `expm1` $e^x - 1$
- `lnp1` $\ln(x+1)$
- `geometric_sum` $\frac{1}{1-x}$
- `sin`
- `cos`
- `tan`
- `sinh`
- `cosh`
- `tanh`
- `asin`
- `atan`
- `asinh`
- `atanh`

Having the capacity of specifying the degree is very important, as users may use other formats than `float64` or `float32` which require higher or lower degree to achieve correct or acceptable precision.

It's possible to define Taylor expansion by implementing a `coeff_at` structure which must meet the following requirement:

- Being template in Integers (`typename`) and index (`size_t`);
- Exposing a type alias `type`, some specialization of `FractionField<Integers>::val`.

For example, to define the serie $1+x+x^2+x^3+\ldots$, users may write:

```cpp
template<typename Integers, size_t i>
struct my_coeff_at {
    using type = typename aerobus::FractionField<Integers>::one;
};

template<typename Integers, size_t degree>
using my_serie = aerobus::taylor<Integers, my_coeff_at, degree>;

static constexpr double x = my_serie<aerobus::i64, 3>::eval(3.0);  // 1 + 3 + 9 + 27
```

On x86-64 and CUDA, the evaluation compiles to a straight chain of multiply-adds, without loop or branch. Fused multiply-adds appear only when the target has them and the compiler may contract `x * y + z`: for example `-mfma` (or `-march=native`) with GCC, which contracts by default with `-std=gnu++20` but not with `-std=c++20` (add `-ffp-contract=fast`). For example, this code:

```cpp
inline void compute_expm1(const size_t N, double* in, double* out) {
    using V = aerobus::expm1<aerobus::i64, 13>;
    for (size_t i = 0; i < N; ++i) {
        out[i] = V::eval(in[i]);
    }
}
```

Yields this assembly (clang 17, `-mavx2 -mfma -O3`) where we can see a pile of Fused-Multiply-Add vector instructions, generated because we unrolled completely the Horner evaluation loop:

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

Given a set (type) satisfies the `IsEuclideanDomain` concept and has no zero divisors, `Aerobus` allows to define its [field of fractions](https://en.wikipedia.org/wiki/Field_of_fractions).

This new type is again a euclidean domain, especially a field, and therefore we can define polynomials over it. The field of fractions of a field is the field itself.

For example, integers modulo a prime `p` form a field, and we can define rational functions over it:

```cpp
using ZpZ = aerobus::zpz<7>;
using Fzpz = aerobus::FractionField<aerobus::polynomial<ZpZ>>;  // rational functions over Z/7Z
```

Integers modulo a composite number have zero divisors (2 * 4 = 0 modulo 8): `FractionField<zpz<8>>` stops with a `static_assert`.

For example, we can easily define [rational functions](https://en.wikipedia.org/wiki/Rational_function) by taking the field of fractions of polynomials:

```cpp
using RF64 = aerobus::FractionField<aerobus::polynomial<aerobus::q64>>;
```

Which also have an evaluation function, as polynomial do.

### Quotient

Given a ring `R`, `Aerobus` provides automatic implementation for [quotient ring](https://en.wikipedia.org/wiki/Quotient_ring) $R/X$ where X is a principal ideal generated by some element, as we know this kind of ideal is two-sided as long as `R` is commutative (and we assume it is).

For example, if we want `R` to be $\mathbb{Z}$ represented as `aerobus::i64`, we can express arithmetic modulo 17 using:

```cpp
using ZpZ = aerobus::Quotient<aerobus::i64, aerobus::i64::val<17>>;
static_assert(ZpZ::is_field);
static_assert(std::is_same_v<ZpZ::inject_constant_t<-1>, ZpZ::inject_constant_t<16>>);
```

As we could have using `zpz<17>`.

Each class has a single type, its canonical representative being the remainder modulo `X` (in `[0, |X|)` for integers). Division multiplies by the inverse modulo `X`. Quotients of integers know whether they are fields; for polynomial rings, pass `true` as third parameter when `X` is irreducible.

This is mainly used to define finite fields of order $p^n$ using Conway polynomials but may have other applications.

## Misc

### Continued Fractions

`Aerobus` gives an implementation for [continued fractions](https://en.wikipedia.org/wiki/Continued_fraction). It can be used this way:

```cpp
using T = aerobus::ContinuedFraction<1, 2, 3, 4>;
constexpr double x = T::val;
constexpr float y = T::get<float>();
```

As practical examples, `aerobus` gives continued fractions of $\pi$, $e$, $\sqrt{2}$ and $\sqrt{3}$:

```cpp
constexpr double A_SQRT3 = aerobus::SQRT3_fraction::val;  // 1.7320508075688772935
```

The value is computed exactly as a `q64` fraction, so long expansions overflow 64 bits (a compile error).

## CUDA

When compiled with `nvcc` and the flag `WITH_CUDA_FP16`, `Aerobus` provides some support of 16 bits floats (`__half` and `__half2`). Compile with `--expt-relaxed-constexpr`.

Unfortunately, NVIDIA did not put enough constexpr in its `cuda_fp16.h` header, so we had to implement our own constexpr static_cast from int16_t to `__half` to make integers polynomials work with `__half`. See [this bug](https://developer.nvidia.com/bugs/4863696).

More, it was not easily possible to make it work for `__half2` because of [another bug](https://developer.nvidia.com/bugs/4872028); with nvcc 13.2, `examples/fp16.cu` compiles without patching `cuda_fp16.h`.

Then nvcc generates splendid assembly, same as for `double` or `float`:

```asm
HFMA2.MMA R5, R6, RZ, 0.0013885498046875, 0.0013885498046875 ;    
HFMA2 R5, R6, R5, 0.008331298828125, 0.008331298828125 ;          
HFMA2.MMA R5, R6, R5, 0.041656494140625, 0.041656494140625 ;      
HFMA2 R5, R6, R5, 0.1666259765625, 0.1666259765625 ;              
HFMA2.MMA R5, R6, R5, 0.5, 0.5 ;                                  
HFMA2 R5, R6, R5, 1, 1 ;                                          
HFMA2.MMA R5, R6, R5, RZ ;                                        
HFMA2 R7, R5, RZ.H0_H0, 0.0013885498046875, 0.0013885498046875 ;  
HFMA2.MMA R7, R5, R7, 0.008331298828125, 0.008331298828125 ;      
HFMA2 R7, R5, R7, 0.041656494140625, 0.041656494140625 ;          
HFMA2.MMA R7, R5, R7, 0.1666259765625, 0.1666259765625 ;          
HFMA2 R7, R5, R7, 0.5, 0.5 ;                                      
HFMA2.MMA R7, R5, R7, 1, 1 ;                                      
HFMA2 R7, R5, R7, RZ.H0_H0 ;                                  
```

Please push to make these bug fixed by NVIDIA.

[![DOI](https://zenodo.org/badge/499577459.svg)](https://zenodo.org/badge/latestdoi/499577459)
