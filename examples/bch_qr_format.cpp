// The BCH code protecting the format information of QR codes, derived at compile time by factoring x^15 - 1.
//
// A QR code stores its error correction level and mask pattern (5 bits) in 15 bits: a BCH(15, 5) codeword
// XORed with 0x5412. The generator of that code is the product of the minimal polynomials of alpha, alpha^3
// and alpha^5 over F2, alpha being a primitive element of GF(16). The minimal polynomials are the irreducible
// factors of x^15 - 1 over F2, which aerobus computes; the compiler then picks those vanishing at alpha^1 ..
// alpha^6 (designed distance 7), checks that their product is the generator of the QR standard
// (x^10 + x^8 + x^5 + x^4 + x^2 + x + 1, 0x537), and builds the 32 format words.
//
// Builds with GCC and clang. Not built with MSVC, like the factorization tests.
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>

#define AEROBUS_CONWAY_IMPORTS
#include "../src/aerobus.h"

// ---------------------------------------------------------------- compile time: the generator

using F2 = aerobus::zpz<2>;
using PF2 = aerobus::polynomial<F2>;

// x^15 - 1, whose roots are the 15 nonzero elements of GF(16)
using x15_minus_1 = PF2::sub_t<PF2::monomial_t<F2::one, 15>, PF2::one>;
using factors = aerobus::factorization_t<x15_minus_1>::factors;

// GF(16) = F2[x] / (x^4 + x + 1), the Conway polynomial CP(2, 4); alpha, the class of x, is primitive
using GF16 = aerobus::GaloisField<2, 4>;
static_assert(std::is_same_v<aerobus::ConwayPolynomial<2, 4>::type, aerobus::make_int_polynomial_t<F2, 1, 0, 0, 1, 1>>);
using alpha = GF16::inject_ring_t<PF2::X>;

// f(alpha^k), for f with coefficients in F2
template<typename f, size_t k, size_t... i>
constexpr bool vanishes_at(std::index_sequence<i...>) {
    using value = aerobus::vadd_t<GF16::mul_t<
        GF16::inject_ring_t<PF2::inject_ring_t<typename f::template coeff_at_t<i>>>, aerobus::pow_t<alpha, i * k>>...>;
    return value::is_zero_v;
}

// f is the minimal polynomial of one of alpha^1, ..., alpha^6
template<typename f, size_t... k>
constexpr bool has_root_among(std::index_sequence<k...>) {
    return (vanishes_at<f, k + 1>(std::make_index_sequence<f::degree + 1>{}) || ...);
}

template<typename List>
struct bch_generator;

template<typename... Fs>
struct bch_generator<aerobus::type_list<Fs...>> {
    using type = aerobus::vmul_t<PF2::one, std::conditional_t<
        has_root_among<typename Fs::type>(std::make_index_sequence<6>{}), typename Fs::type, PF2::one>...>;
};

using generator = bch_generator<factors>::type;

// x^10 + x^8 + x^5 + x^4 + x^2 + x + 1, the generator in ISO/IEC 18004
static_assert(std::is_same_v<generator, aerobus::make_int_polynomial_t<F2, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1>>,
              "the product of the minimal polynomials of alpha .. alpha^6 is the QR format generator");

// ---------------------------------------------------------------- compile time: the format words

// a polynomial over F2 as bits: bit i is the coefficient of x^i
template<typename P, size_t... i>
constexpr uint32_t to_bits(std::index_sequence<i...>) {
    return ((static_cast<uint32_t>(P::template coeff_at_t<i>::v) << i) | ...);
}

constexpr uint32_t G = to_bits<generator>(std::make_index_sequence<generator::degree + 1>{});
static_assert(G == 0x537);

constexpr uint32_t FORMAT_MASK = 0x5412;

// systematic encoding: data x^10 + (data x^10 mod g), then the mask
constexpr uint16_t encode_format(uint32_t data) {
    uint32_t remainder = data << 10;
    for (int bit = 14; bit >= 10; --bit) {
        if ((remainder >> bit) & 1U) {
            remainder ^= G << (bit - 10);
        }
    }
    return static_cast<uint16_t>(((data << 10) | remainder) ^ FORMAT_MASK);
}

constexpr std::array<uint16_t, 32> format_words() {
    std::array<uint16_t, 32> words{};
    for (uint32_t data = 0; data < 32; ++data) {
        words[data] = encode_format(data);
    }
    return words;
}

constexpr auto FORMAT = format_words();

// the table of ISO/IEC 18004, indexed by (error correction level bits << 3) | mask, with L = 01, M = 00,
// Q = 11, H = 10
constexpr std::array<uint16_t, 32> STANDARD = {
    0x5412, 0x5125, 0x5e7c, 0x5b4b, 0x45f9, 0x40ce, 0x4f97, 0x4aa0,  // M, masks 0 to 7
    0x77c4, 0x72f3, 0x7daa, 0x789d, 0x662f, 0x6318, 0x6c41, 0x6976,  // L
    0x1689, 0x13be, 0x1ce7, 0x19d0, 0x0762, 0x0255, 0x0d0c, 0x083b,  // H
    0x355f, 0x3068, 0x3f31, 0x3a06, 0x24b4, 0x2183, 0x2eda, 0x2bed,  // Q
};
static_assert(FORMAT == STANDARD, "the 32 format words are those of the QR standard");

constexpr int popcount15(uint32_t v) {
    int n = 0;
    for (int bit = 0; bit < 15; ++bit) {
        n += static_cast<int>((v >> bit) & 1U);
    }
    return n;
}

constexpr int minimum_distance() {
    int best = 15;
    for (size_t i = 0; i < 32; ++i) {
        for (size_t j = i + 1; j < 32; ++j) {
            const int d = popcount15(static_cast<uint32_t>(FORMAT[i] ^ FORMAT[j]));
            best = d < best ? d : best;
        }
    }
    return best;
}
static_assert(minimum_distance() == 7, "BCH bound: designed distance 7, so 3 errors are corrected");

// ---------------------------------------------------------------- run time: decoding

// the data bits of the nearest format word, as a QR reader does
static uint32_t decode_format(uint16_t received) {
    uint32_t best = 0;
    int best_distance = 16;
    for (uint32_t data = 0; data < 32; ++data) {
        const int d = popcount15(static_cast<uint32_t>(received ^ FORMAT[data]));
        if (d < best_distance) {
            best_distance = d;
            best = data;
        }
    }
    return best;
}

static int failures = 0;

static void check(bool condition, const char* what) {
    std::printf("%s %s\n", condition ? "ok  " : "FAIL", what);
    failures += condition ? 0 : 1;
}

template<typename... Fs>
static void print_factors(aerobus::type_list<Fs...>) {
    ((std::printf("  %s%s\n", Fs::type::to_string().c_str(),
                  has_root_among<typename Fs::type>(std::make_index_sequence<6>{}) ? "   <- in the generator" : "")),
     ...);
}

int main() {
    std::printf("x^15 - 1 over F2 factors as:\n");
    print_factors(factors{});
    std::printf("generator: %s\n", generator::to_string().c_str());

    check(true, "generator is x^10 + x^8 + x^5 + x^4 + x^2 + x + 1 (static_assert)");
    check(true, "the 32 format words match the QR standard (static_assert)");
    check(true, "minimum distance is 7 (static_assert)");

    // every word, with every pattern of up to 3 flipped bits among the 15
    int corrected = 0;
    int patterns = 0;
    for (uint32_t data = 0; data < 32; ++data) {
        for (uint32_t error = 0; error < (1U << 15); ++error) {
            if (popcount15(error) > 3) {
                continue;
            }
            ++patterns;
            corrected += decode_format(static_cast<uint16_t>(FORMAT[data] ^ error)) == data ? 1 : 0;
        }
    }
    std::printf("     %d of %d corrupted format words decoded\n", corrected, patterns);
    check(corrected == patterns && patterns == 32 * 576, "every format word with up to 3 flipped bits is decoded");

    return failures == 0 ? 0 : 1;
}
