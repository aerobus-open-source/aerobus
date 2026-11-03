# Reed-Solomon with a field the compiler builds and proves

November 2026

Every QR code carries Reed-Solomon error correction over GF(256), the field with 256 elements. This post builds that field with [aerobus](https://github.com/aerobus-open-source/aerobus), a header-only C++20 library where every algebraic value is a type: the compiler computes the exp/log tables and the generator polynomial, and a `static_assert` proves the table is a generator of the field before the program exists. At run time, a codec of about 110 lines only reads constants, and it reproduces the error correction bytes of the QR standard's "HELLO WORLD" example.

## The field is a type

GF(256) is the ring of polynomials over Z/2Z, taken modulo an irreducible polynomial of degree 8. In aerobus that sentence is one line:

```cpp
using GF256 = aerobus::GaloisField<2, 8>;  // F2[x] / CP(2, 8)
```

`GaloisField<p, n>` is `Quotient<polynomial<zpz<p>>, ConwayPolynomial<p, n>::type, true>`. Values are types: `GF256::mul_t<a, b>` is the product, computed by the compiler as polynomial multiplication followed by a remainder.

The modulus is the Conway polynomial CP(2, 8) = x⁸ + x⁴ + x³ + x² + 1, which reads 0x11D in binary. That is exactly the polynomial of the QR code standard (and of Data Matrix and RAID-6). Nothing is tuned to match QR: the canonical choice of number theory and the choice of the standard coincide, so the tables below come out right by construction.

## Tables computed, then proven, by the compiler

Reed-Solomon arithmetic runs on two tables: exp[i] = αⁱ and its inverse, log. Here α is the class of x. The compiler computes its powers by repeated multiplication in the field:

```cpp
using alpha = GF256::inject_ring_t<aerobus::polynomial<aerobus::zpz<2>>::X>;

template<size_t i>
struct alpha_power {
    using type = GF256::mul_t<typename alpha_power<i - 1>::type, alpha>;
};

template<>
struct alpha_power<0> {
    using type = GF256::one;
};
```

Each element is a polynomial of degree below 8 over Z/2Z; its coefficients are the bits of a byte:

```cpp
template<typename E, size_t... bit>
constexpr uint8_t pack(std::index_sequence<bit...>) {
    using P = typename E::type;  // canonical representative, a polynomial
    return static_cast<uint8_t>(((P::template coeff_at_t<bit>::v << bit) | ...));
}

template<size_t... i>
constexpr std::array<uint8_t, 255> exp_table(std::index_sequence<i...>) {
    return {pack<typename alpha_power<i>::type>(std::make_index_sequence<8>{})...};
}

constexpr auto EXP = exp_table(std::make_index_sequence<255>{});
```

The log table, the inverses and the full 256 × 256 product table are short `constexpr` loops over EXP. Then the compiler checks the property everything else relies on: α generates the 255 nonzero bytes.

```cpp
constexpr bool alpha_is_primitive() {
    std::array<bool, 256> seen{};
    for (uint8_t value : EXP) {
        if (value == 0 || seen[value]) return false;
        seen[value] = true;
    }
    return true;
}
static_assert(alpha_is_primitive(), "alpha must generate the multiplicative group of GF(256)");
static_assert(EXP[8] == 0x1D, "alpha^8 = x^4 + x^3 + x^2 + 1");
```

If the modulus were not primitive, or a table entry wrong, the program would not compile. The tables land in `.rodata`: no initialisation at start-up, no static constructor.

A plain `constexpr` loop shifting bytes would also build the table. The difference is where the definition lives: here the code *is* the mathematical definition (powers of a generator in F₂[x]/(CP(2, 8))), and the same two lines give F₉, F₂₅ or any F_pⁿ from the Conway table.

## The generator polynomial is a type too

A Reed-Solomon code with n parity bytes is defined by its generator g(x) = (x − α⁰)(x − α¹)…(x − αⁿ⁻¹). Polynomials over GF(256) are `polynomial<GF256>`, so g is a product of types (−αⁱ = αⁱ in characteristic 2):

```cpp
using PGF = aerobus::polynomial<GF256>;

template<size_t n>
struct rs_generator {
    using type = PGF::mul_t<typename rs_generator<n - 1>::type,
                            PGF::val<GF256::one, typename alpha_power<n - 1>::type>>;
};

template<>
struct rs_generator<0> {
    using type = PGF::one;
};
```

For the 10 parity bytes of a version 1-M QR code, the compiler finds 1, 216, 194, 159, 111, 199, 94, 95, 113, 157, 193, which is α⁰, α²⁵¹, α⁶⁷, …, α⁴⁵, the generator printed in QR code references.

## Encoding and decoding "HELLO WORLD"

At run time there is no aerobus at all: a systematic encoder (long division by g) and a textbook decoder (syndromes, Berlekamp-Massey, Chien search, Forney) read the tables and the generator bytes. The example checks them against an independent implementation (the Python library `reedsolo`) and the QR standard:

| Check | Result |
| --- | --- |
| exp table starts 01 02 04 … 80 1d 3a | ok |
| α generates the 255 nonzero bytes (also a `static_assert`) | ok |
| product table equals carry-less multiplication for the 65,536 pairs | ok |
| generator for 10 parity bytes | ok |
| "HELLO WORLD" data codewords 32 91 11 120 … encode to 196 35 39 119 235 215 231 226 93 23 | ok |
| 5 corrupted bytes in that block are corrected | ok |
| 10,000 random codewords with up to 5 corrupted bytes | 10,000 corrected |

## How fast is it?

With its multiplication table also computed at compile time, the example codec is the fastest of the three everywhere except decoding long blocks, where Schifra stays ahead. The speed belongs to the codec, not to aerobus: aerobus contributes the tables, their proof, and zero start-up cost.

The comparison uses two reference implementations: the Reed-Solomon code of [zxing-cpp](https://github.com/zxing-cpp/zxing-cpp) v2.3.0, the C++ barcode library, and [Schifra](https://github.com/ArashPartow/schifra), a dedicated Reed-Solomon library. On 1,000 random blocks per case, all three produce the same parity bytes and correct the same errors. Times are the median of 7 repeats; two runs agreed within a few percent. AMD Ryzen 5 PRO 5675U, GCC 15.2, `-O2`.

**QR version 1-M block: 16 data + 10 parity bytes, 5 errors** (ns per block)

| Operation | aerobus example | zxing-cpp | Schifra |
| --- | --- | --- | --- |
| Encode | 161 | 300 (960 through `ReedSolomonEncode`) | 3,950 |
| Decode, 5 errors | 965 | 3,160 | 6,000 to 6,400 |
| Decode, no error | 210 | 670 | 2,810 |

**RS(255,223): 223 data + 32 parity bytes, 16 errors** (ns per block)

| Operation | aerobus example | zxing-cpp | Schifra |
| --- | --- | --- | --- |
| Encode | 3,500 | 7,300 | 11,300 |
| Decode, 16 errors | 21,000 | 40,900 | 16,400 to 16,900 |
| Decode, no error | 10,700 | 28,700 | 6,500 |

**Table set-up, once per process:** aerobus 0 ns (constants), zxing-cpp 3.7 µs (first use of its static field), Schifra 426 µs (field and generator).

What the numbers say:

- The full 256 × 256 product table (64 KiB) is `constexpr` data derived from the proven exp/log tables, one constant per row so that each compile-time evaluation stays within clang's and MSVC's step limits. It halved or better every time of the first version of the example, which multiplied through log/exp with a modulo.
- On short blocks the example is 2 to 3 times faster than zxing-cpp and 6 to 25 times faster than Schifra, which pads short codes to 255 bytes. zxing-cpp's convenience function `ReedSolomonEncode` rebuilds the generator at every call.
- On long blocks Schifra still decodes faster (about 1.3 times with errors, 1.6 times without). It pays 430 µs at start-up for its tables, where the example pays nothing at run time and about 4 s at compile time.

## Limits

- **Compile time.** The example compiles in about 4 s with GCC and 7 s with clang; the 254 field products, the degree-10 generator and the 64 KiB product table are the cost. A generator with 32 parity bytes takes longer.
- **Not with MSVC.** GCC and clang compile the example; MSVC did not finish within 13 minutes in CI, so the example is not built with MSVC.
- **Template depth.** The recursion over 255 powers needs `-ftemplate-depth=30000` with GCC and clang, as set by the example's CMake target.
- **Not for cryptography.** Table lookups are not constant-time. Use this for error correction, not for AES or secret sharing in production.
- **Speed is the decoder's.** Aerobus builds and checks the field; the run-time throughput is whatever decoder you write on top. Schifra's decoder is still faster on 255-byte blocks.
- **Same tables as a loop.** Ten lines of `constexpr` bit shifting also produce `EXP`. The point here is that the definition, the proof and every other field come from the same few lines.

## Try it

The exp table, its proof and the generator run in your browser on [Compiler Explorer](https://godbolt.org/clientstate/eyJzZXNzaW9ucyI6W3siaWQiOjEsImxhbmd1YWdlIjoiYysrIiwic291cmNlIjoiLy8gR0YoMjU2KSwgdGhlIGZpZWxkIG9mIFFSIGNvZGVzLCBidWlsdCBhbmQgcHJvdmVuIGJ5IHRoZSBjb21waWxlciB3aXRoIGFlcm9idXMuXG4vLyBBcnRpY2xlOiBodHRwczovL2Flcm9idXMtb3Blbi1zb3VyY2UuZ2l0aHViLmlvL2Flcm9idXMvbWRfYXJ0aWNsZXNfMnJlZWRfX3NvbG9tb25fX3FyLmh0bWxcbiNpbmNsdWRlIDxhcnJheT5cbiNpbmNsdWRlIDxjc3RkaW50PlxuI2luY2x1ZGUgPGNzdGRpbz5cbiNpbmNsdWRlIDx1dGlsaXR5PlxuXG4jZGVmaW5lIEFFUk9CVVNfQ09OV0FZX0lNUE9SVFNcbiNpbmNsdWRlIDxodHRwczovL3Jhdy5naXRodWJ1c2VyY29udGVudC5jb20vYWVyb2J1cy1vcGVuLXNvdXJjZS9hZXJvYnVzLzM0OGEzM2YzNDIyZjM0ZTc3MjgzMTEyMGYwNDEzZWQ3MmEzMmI5OTcvc3JjL2Flcm9idXMuaD5cblxudXNpbmcgR0YyNTYgPSBhZXJvYnVzOjpHYWxvaXNGaWVsZDwyLCA4PjsgIC8vIEYyW3hdIC8gKHheOCArIHheNCArIHheMyArIHheMiArIDEpXG51c2luZyBhbHBoYSA9IEdGMjU2OjppbmplY3RfcmluZ190PGFlcm9idXM6OnBvbHlub21pYWw8YWVyb2J1czo6enB6PDI%2BPjo6WD47XG5cbi8vIGFscGhhXmksIGNvbXB1dGVkIGJ5IHRoZSBjb21waWxlclxudGVtcGxhdGU8c2l6ZV90IGk%2BXG5zdHJ1Y3QgYWxwaGFfcG93ZXIge1xuICAgIHVzaW5nIHR5cGUgPSBHRjI1Njo6bXVsX3Q8dHlwZW5hbWUgYWxwaGFfcG93ZXI8aSAtIDE%2BOjp0eXBlLCBhbHBoYT47XG59O1xuXG50ZW1wbGF0ZTw%2BXG5zdHJ1Y3QgYWxwaGFfcG93ZXI8MD4ge1xuICAgIHVzaW5nIHR5cGUgPSBHRjI1Njo6b25lO1xufTtcblxuLy8gYW4gZWxlbWVudCBhcyBhIGJ5dGU6IGJpdCBpIGlzIHRoZSBjb2VmZmljaWVudCBvZiB4XmlcbnRlbXBsYXRlPHR5cGVuYW1lIEUsIHNpemVfdC4uLiBiaXQ%2BXG5jb25zdGV4cHIgdWludDhfdCBwYWNrKHN0ZDo6aW5kZXhfc2VxdWVuY2U8Yml0Li4uPikge1xuICAgIHVzaW5nIFAgPSB0eXBlbmFtZSBFOjp0eXBlO1xuICAgIHJldHVybiBzdGF0aWNfY2FzdDx1aW50OF90PigoKFA6OnRlbXBsYXRlIGNvZWZmX2F0X3Q8Yml0Pjo6diA8PCBiaXQpIHwgLi4uKSk7XG59XG5cbnRlbXBsYXRlPHNpemVfdC4uLiBpPlxuY29uc3RleHByIHN0ZDo6YXJyYXk8dWludDhfdCwgMjU1PiBleHBfdGFibGUoc3RkOjppbmRleF9zZXF1ZW5jZTxpLi4uPikge1xuICAgIHJldHVybiB7cGFjazx0eXBlbmFtZSBhbHBoYV9wb3dlcjxpPjo6dHlwZT4oc3RkOjptYWtlX2luZGV4X3NlcXVlbmNlPDg%2Be30pLi4ufTtcbn1cblxuY29uc3RleHByIGF1dG8gRVhQID0gZXhwX3RhYmxlKHN0ZDo6bWFrZV9pbmRleF9zZXF1ZW5jZTwyNTU%2Be30pO1xuXG4vLyBwcm92ZW4gYXQgY29tcGlsZSB0aW1lOiBhbHBoYSBnZW5lcmF0ZXMgdGhlIDI1NSBub256ZXJvIGJ5dGVzXG5jb25zdGV4cHIgYm9vbCBhbHBoYV9pc19wcmltaXRpdmUoKSB7XG4gICAgc3RkOjphcnJheTxib29sLCAyNTY%2BIHNlZW57fTtcbiAgICBmb3IgKHVpbnQ4X3QgdmFsdWUgOiBFWFApIHtcbiAgICAgICAgaWYgKHZhbHVlID09IDAgfHwgc2Vlblt2YWx1ZV0pIHJldHVybiBmYWxzZTtcbiAgICAgICAgc2Vlblt2YWx1ZV0gPSB0cnVlO1xuICAgIH1cbiAgICByZXR1cm4gdHJ1ZTtcbn1cbnN0YXRpY19hc3NlcnQoYWxwaGFfaXNfcHJpbWl0aXZlKCksIFwiYWxwaGEgbXVzdCBnZW5lcmF0ZSB0aGUgbXVsdGlwbGljYXRpdmUgZ3JvdXAgb2YgR0YoMjU2KVwiKTtcblxuLy8gdGhlIFJlZWQtU29sb21vbiBnZW5lcmF0b3IgKHggLSBhXjApLi4uKHggLSBhXjkpIG9mIGEgdmVyc2lvbiAxLU0gUVIgY29kZSwgYXMgYSBwb2x5bm9taWFsIG92ZXIgR0YoMjU2KVxudXNpbmcgUEdGID0gYWVyb2J1czo6cG9seW5vbWlhbDxHRjI1Nj47XG5cbnRlbXBsYXRlPHNpemVfdCBuPlxuc3RydWN0IHJzX2dlbmVyYXRvciB7XG4gICAgdXNpbmcgdHlwZSA9IFBHRjo6bXVsX3Q8dHlwZW5hbWUgcnNfZ2VuZXJhdG9yPG4gLSAxPjo6dHlwZSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBQR0Y6OnZhbDxHRjI1Njo6b25lLCB0eXBlbmFtZSBhbHBoYV9wb3dlcjxuIC0gMT46OnR5cGU%2BPjtcbn07XG5cbnRlbXBsYXRlPD5cbnN0cnVjdCByc19nZW5lcmF0b3I8MD4ge1xuICAgIHVzaW5nIHR5cGUgPSBQR0Y6Om9uZTtcbn07XG5cbnRlbXBsYXRlPHR5cGVuYW1lIEcsIHNpemVfdC4uLiBpPlxuY29uc3RleHByIHN0ZDo6YXJyYXk8dWludDhfdCwgc2l6ZW9mLi4uKGkpPiBjb2VmZmljaWVudHMoc3RkOjppbmRleF9zZXF1ZW5jZTxpLi4uPikge1xuICAgIHJldHVybiB7cGFjazx0eXBlbmFtZSBHOjp0ZW1wbGF0ZSBjb2VmZl9hdF90PEc6OmRlZ3JlZSAtIGk%2BPihzdGQ6Om1ha2VfaW5kZXhfc2VxdWVuY2U8OD57fSkuLi59O1xufVxuXG5jb25zdGV4cHIgYXV0byBHMTAgPSBjb2VmZmljaWVudHM8cnNfZ2VuZXJhdG9yPDEwPjo6dHlwZT4oc3RkOjptYWtlX2luZGV4X3NlcXVlbmNlPDExPnt9KTtcbnN0YXRpY19hc3NlcnQoRzEwID09IHN0ZDo6YXJyYXk8dWludDhfdCwgMTE%2BezEsIDIxNiwgMTk0LCAxNTksIDExMSwgMTk5LCA5NCwgOTUsIDExMywgMTU3LCAxOTN9KTtcblxuaW50IG1haW4oKSB7XG4gICAgc3RkOjpwcmludGYoXCJleHAgdGFibGU6XCIpO1xuICAgIGZvciAoaW50IGkgPSAwOyBpIDwgMTY7ICsraSkgc3RkOjpwcmludGYoXCIgJTAyeFwiLCBFWFBbaV0pO1xuICAgIHN0ZDo6cHJpbnRmKFwiIC4uLlxcbmdlbmVyYXRvciBmb3IgMTAgcGFyaXR5IGJ5dGVzOlwiKTtcbiAgICBmb3IgKHVpbnQ4X3QgYyA6IEcxMCkgc3RkOjpwcmludGYoXCIgJWRcIiwgYyk7XG4gICAgc3RkOjpwcmludGYoXCJcXG5cIik7XG59XG4iLCJjb21waWxlcnMiOlt7ImlkIjoiZzE1MiIsIm9wdGlvbnMiOiItc3RkPWMrKzIwIC1PMiAtZnRlbXBsYXRlLWRlcHRoPTMwMDAwIn1dLCJleGVjdXRvcnMiOlt7ImNvbXBpbGVyIjp7ImlkIjoiZzE1MiIsIm9wdGlvbnMiOiItc3RkPWMrKzIwIC1PMiAtZnRlbXBsYXRlLWRlcHRoPTMwMDAwIiwibGlicyI6W119LCJhcmd1bWVudHMiOiIiLCJzdGRpbiI6IiJ9XX1dfQ%3D%3D) (GCC 15.2; the header is included from GitHub at the commit of this article).

The whole example is [`examples/reed_solomon.cpp`](https://github.com/aerobus-open-source/aerobus/blob/main/examples/reed_solomon.cpp) in the aerobus repository. It runs as one of the library's tests:

```bash
git clone https://github.com/aerobus-open-source/aerobus.git
cd aerobus
cmake -S . -B build
cmake --build build --target example_reed_solomon
./build/examples/example_reed_solomon
```

Aerobus is a single header (`src/aerobus.h`, C++20, GCC 13+, clang 17+, MSVC); define `AEROBUS_CONWAY_IMPORTS` to get the Conway polynomials, hence `GaloisField<p, n>` for every prime p below 1000. The code blocks of this article are compiled by the repository's tests.
