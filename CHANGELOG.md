# Changelog

## Unreleased

### Added

- `bigint`: arbitrary precision integers as an euclidean domain, with `from_string_t`, `from_hex_t`, `to_string`, `to_hex` and `get`. The arithmetic runs in constexpr functions, division is Knuth's algorithm D and gcd a constexpr Euclid. Predefined families and Taylor series accept `bigint` integers, without degree limit. `string_literal` for string template arguments.
- Finite fields and factorization: concept `IsFiniteField`, `field_order_v`, `characteristic_v`, `field_element_t`, `monic_t`, `is_irreducible_v` (Rabin), `squarefree_t` over any field, and `factorization_t` over finite fields (distinct degree factorization and a deterministic Cantor-Zassenhaus), returning a `factorization` of sorted monic `factor`s.
- Example `examples/bch_qr_format.cpp`: the BCH(15, 5) generator of the QR format information, derived by factoring x^15 - 1 over F2, and the 32 format words of the standard. Not built with MSVC.
- Article [The QR format code is a factor of x¹⁵ − 1](articles/bch_qr_format.md), whose code blocks and quoted numbers are compiled by CTest.
- Concept `coefficients_fit<Family, deg, I>`, required by every predefined polynomial family and Taylor series. A degree whose coefficients overflow `i16`, `i32` or `i64` is now a constraint failure naming the concept, instead of an overflow error inside the header, and can be tested with a requires-expression. The limits were measured; `tests/degree_limits_tests.cpp` checks that those of the Taylor series are exact.

### Changed

- `factorial`, `combination`, `bernoulli` and `alternate` work for any integer ring (`inject_constant_t` instead of `val<x>`); the `value` members of their internal structs are gone, use `factorial_v`, `combination_v` and `alternate_v`.

## 2.0.0 (2026-11-03)

Version 2.0 fixes the defects found by the code review of commit `f7e0c79` (November 2026). It breaks the API of 1.x.

### Breaking changes

| 1.x | 2.0 |
| --- | --- |
| `aerobus::div_t<X, Y>` (ambiguous with the C type `div_t`) | `aerobus::divide_t<X, Y>` |
| `INLINED`, `DEVICE`, `ALIGNED` macros | `AEROBUS_INLINE`, `AEROBUS_DEVICE`, `AEROBUS_ALIGNED` |
| `ZPZV`, `POLYV` defined after the include | undefined after the Conway table |
| `pow_t<T, p, n>`, `pow_v<T, p, n>` | `pow_t<p, n>`, `pow_v<p, n>` (ring inferred from `p`) |
| `gcd_t<T, A, B>` | `gcd_t<A, B>` (ring inferred from `A`) |
| `type_list::insert<T, index>` | `type_list::insert<index, T>` |
| `int` indices of `alternate_t`, `stirling_*` | `size_t` indices |
| `Quotient<Ring, X>::val<V>::raw_t` (unreduced) | `type`, the canonical representative |
| `Quotient<Ring, X>` requires `IsRing` | requires `IsEuclideanDomain`; third parameter `declared_field` |
| `polynomial<Ring>::horner_reduction_t` public | private (use `val::value_at_t`) |
| `struct i16`, `i32`, `i64` | aliases of `integers<int16_t>`, `integers<int32_t>`, `integers<int64_t>` |
| `i16`, `q16` only with `WITH_CUDA_FP16` | always defined |
| `FractionField<R>` for any euclidean `R` | `R` must have no zero divisors (`zpz<8>` is rejected) |
| `polynomial<Ring>::val<...>` kept leading zeros | removes them: `val<0, 1, 2>` is `val<1, 2>` |
| `polynomial<Ring>::val` accepted coefficients of any ring | static_asserts they are values of `Ring` |
| `zpz<p>::val<x>` stored `x % p` (possibly negative) | alias of `residue<x mod p>` in `[0, p)` |
| `zpz<p>::div_t`, `mod_t`: integer division of representatives | multiplication by the inverse; remainder 0 |
| `_FractionField` (reserved name) | `internal::fraction_field_impl` |
| `aerobus::libm` leftovers in `src/` and `examples/` | moved to `attic/` |

After `using namespace aerobus`, the Taylor series `exp`, `sin`, `cos`, `tan`, `atanh`... are ambiguous with the C functions of `<math.h>`: qualify them with `aerobus::`.

### Fixed

- `is_prime` reported 35, 49, 77, 143 as prime and did not compile for 25, 55, 121, 187.
- `zpz<p>`: negative representatives, wrong comparisons, `int32_t` overflow for large `p`, division by representatives, `gcd_t` returning an `i32` value.
- `Quotient`: the canonical form negated negative remainders (`-1 + 1 != 0` in Z/17Z, wrong products in F9); division divided representatives; `sub_t`, `gcd_t`, `gt_t`, `lt_t`, `is_field` were missing; `Embed` returned the unreduced value.
- `IsEuclideanDomain` and `IsField` did not check their flags (`IsField<i32>` was true).
- Integer evaluation used swapped `int16_t`/`int32_t` helpers (`x^2` at 1000 gave 16960).
- Compensated Horner: the error-free product was wrong for `float` (split constant) and on clang (no fused multiply-add); it now uses `std::fma`.
- Horner evaluation returned NaN at infinity.
- Fraction fields: order relations (`gt_t`) wrong or not compiling, `pos_t`, `gcd_t`, division by zero, rational functions equal up to a constant had different canonical forms.
- Polynomials: zero tests by type identity, products not simplified over rings with zero divisors, endless recursion on inexact division, `gt_t` false for equal degrees, zero printed as an empty string.
- `asin`, `tan`, `tanh` Taylor series did not compile.
- `combination` and `stirling_2` with `k > n`; `touchard` mixing rings; `hermite_phys` building `pi64` for any ring; `i64::val::inner_type` was `int32_t`.
- `aligned_malloc` without `aligned_free` and with an invalid size for `aligned_alloc`.
- CUDA: `arithmetic_helpers<__half>::two_pi` was negative; `meta_libm<__half>::floor` was device only.

### Added

- `GaloisField<p, n>`, `Quotient<Ring, X, declared_field>`, `integers<Int>`, `aligned_free`, `ContinuedFraction::get<T>()`, `AEROBUS_VERSION_MAJOR/MINOR/PATCH`, `is_integral_domain` on every ring, the same member set on every built-in ring.
- CMake: `aerobus::aerobus` target, install rules and `find_package(Aerobus)`, options `AEROBUS_BUILD_TESTS`, `AEROBUS_BUILD_EXAMPLES`, `AEROBUS_BUILD_BENCHMARKS`; OpenMP only for the benchmarks.
- Tests split by topic; compile-fail tests for every static assertion; ring axiom tests; examples, README blocks, the generated Conway table and the CUDA fp16 host code checked by CTest.
- `tools/generate_conway.py`, which regenerates the Conway table.
- Example `examples/reed_solomon.cpp`: a Reed-Solomon codec over `GaloisField<2, 8>` (the QR code field) whose tables and generator are computed and checked at compile time, with the article [Reed-Solomon with a field the compiler builds and proves](articles/reed_solomon_qr.md), whose code blocks are compiled by CTest. Not built with MSVC, which does not finish compiling it.

### Changed

- CI: GCC 14, clang 18, macOS and MSVC test jobs, sanitizers, lint (cpplint, cppcheck, Conway table check) and a documentation job whose Doxygen warnings are errors; pushes are built on `main` only, pull requests once. The paper build is removed.

## 1.5.2 and earlier

See the git tags.
