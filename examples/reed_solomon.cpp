// Reed-Solomon error correction of QR codes, with GF(256) built by aerobus at compile time.
//
// GF(256) is GaloisField<2, 8>: polynomials over Z/2Z modulo the Conway polynomial
// x^8 + x^4 + x^3 + x^2 + 1 (0x11D), the field of QR codes, Data Matrix and RAID-6.
// The compiler computes the powers of alpha (the class of x), hence the exp/log tables, and the
// generator polynomial (x - alpha^0)(x - alpha^1)...(x - alpha^9) as a value of polynomial<GF(256)>.
// At run time, encoding and decoding only read the tables.
#include <array>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

#define AEROBUS_CONWAY_IMPORTS
#include "../src/aerobus.h"

// ---------------------------------------------------------------- compile time: tables

using F2 = aerobus::zpz<2>;
using PF2 = aerobus::polynomial<F2>;
using GF256 = aerobus::GaloisField<2, 8>;
using alpha = GF256::inject_ring_t<PF2::X>;  // the class of x, a generator of the multiplicative group

// alpha^i, computed by the compiler
template<size_t i>
struct alpha_power {
    using type = GF256::mul_t<typename alpha_power<i - 1>::type, alpha>;
};

template<>
struct alpha_power<0> {
    using type = GF256::one;
};

// an element of GF(256) as a byte: bit i is the coefficient of x^i of its representative
template<typename E, size_t... bit>
constexpr uint8_t pack(std::index_sequence<bit...>) {
    using P = typename E::type;
    return static_cast<uint8_t>(((P::template coeff_at_t<bit>::v << bit) | ...));
}

template<typename E>
constexpr uint8_t to_byte() {
    return pack<E>(std::make_index_sequence<8>{});
}

template<size_t... i>
constexpr std::array<uint8_t, 255> exp_table(std::index_sequence<i...>) {
    return {to_byte<typename alpha_power<i>::type>()...};
}

// exp[i] = alpha^i
constexpr std::array<uint8_t, 255> exp_table() {
    return exp_table(std::make_index_sequence<255>{});
}

// log[alpha^i] = i (log[0] is unused)
constexpr std::array<uint8_t, 256> log_table() {
    std::array<uint8_t, 256> log{};
    const auto exp = exp_table();
    for (size_t i = 0; i < exp.size(); ++i) {
        log[exp[i]] = static_cast<uint8_t>(i);
    }
    return log;
}

// generator polynomial of Reed-Solomon codes with n parity bytes: (x - alpha^0)...(x - alpha^(n-1)),
// a value of polynomial<GF(256)>; -alpha^i = alpha^i in characteristic 2
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

template<typename G, size_t... i>
constexpr std::array<uint8_t, sizeof...(i)> coefficient_bytes(std::index_sequence<i...>) {
    return {to_byte<typename G::template coeff_at_t<G::degree - i>>()...};
}

// coefficients of the generator as bytes, highest degree first
template<size_t parity>
constexpr std::array<uint8_t, parity + 1> generator_bytes() {
    using G = typename rs_generator<parity>::type;
    static_assert(G::degree == parity);
    return coefficient_bytes<G>(std::make_index_sequence<parity + 1>{});
}

constexpr auto EXP = exp_table();
constexpr auto LOG = log_table();

// checked by the compiler: alpha has order 255, its powers are the 255 nonzero bytes
constexpr bool alpha_is_primitive() {
    std::array<bool, 256> seen{};
    for (uint8_t value : EXP) {
        if (value == 0 || seen[value]) {
            return false;
        }
        seen[value] = true;
    }
    return true;
}
static_assert(alpha_is_primitive(), "alpha must generate the multiplicative group of GF(256)");
static_assert(EXP[8] == 0x1D, "alpha^8 = x^4 + x^3 + x^2 + 1, the Conway polynomial 0x11D");

// ---------------------------------------------------------------- run time: codec

static uint8_t gf_mul(uint8_t a, uint8_t b) {
    return (a == 0 || b == 0) ? 0 : EXP[(LOG[a] + LOG[b]) % 255];
}

static uint8_t gf_div(uint8_t a, uint8_t b) {
    return a == 0 ? 0 : EXP[(LOG[a] + 255 - LOG[b]) % 255];
}

static uint8_t gf_pow_alpha(int e) {
    return EXP[((e % 255) + 255) % 255];
}

// value at x of a polynomial stored lowest degree first
static uint8_t eval_low_first(const std::vector<uint8_t>& p, uint8_t x) {
    uint8_t y = 0;
    for (size_t i = p.size(); i-- > 0;) {
        y = gf_mul(y, x) ^ p[i];
    }
    return y;
}

// systematic encoding: data followed by the remainder of data * x^parity divided by the generator
template<size_t parity>
std::vector<uint8_t> rs_encode(const std::vector<uint8_t>& data) {
    constexpr auto g = generator_bytes<parity>();
    std::vector<uint8_t> remainder(parity, 0);
    for (uint8_t byte : data) {
        const uint8_t factor = byte ^ remainder[0];
        remainder.erase(remainder.begin());
        remainder.push_back(0);
        for (size_t j = 0; j < parity; ++j) {
            remainder[j] ^= gf_mul(g[j + 1], factor);
        }
    }
    std::vector<uint8_t> codeword = data;
    codeword.insert(codeword.end(), remainder.begin(), remainder.end());
    return codeword;
}

// corrects up to parity / 2 byte errors in place; false when the codeword cannot be corrected
template<size_t parity>
bool rs_decode(std::vector<uint8_t>* codeword) {
    std::vector<uint8_t>& c = *codeword;
    const int n = static_cast<int>(c.size());

    // syndromes S_j = c(alpha^j), c[0] being the coefficient of x^(n-1)
    std::vector<uint8_t> syndromes(parity);
    bool clean = true;
    for (size_t j = 0; j < parity; ++j) {
        uint8_t s = 0;
        for (int k = 0; k < n; ++k) {
            s = gf_mul(s, gf_pow_alpha(static_cast<int>(j))) ^ c[k];
        }
        syndromes[j] = s;
        clean = clean && s == 0;
    }
    if (clean) {
        return true;
    }

    // Berlekamp-Massey: error locator lambda, lowest degree first
    std::vector<uint8_t> lambda = {1}, previous = {1};
    size_t errors = 0;
    int shift = 1;
    uint8_t previous_discrepancy = 1;
    for (size_t step = 0; step < parity; ++step) {
        uint8_t discrepancy = syndromes[step];
        for (size_t i = 1; i <= errors && i < lambda.size(); ++i) {
            discrepancy ^= gf_mul(lambda[i], syndromes[step - i]);
        }
        if (discrepancy == 0) {
            ++shift;
            continue;
        }
        std::vector<uint8_t> updated = lambda;
        const uint8_t scale = gf_div(discrepancy, previous_discrepancy);
        if (updated.size() < previous.size() + shift) {
            updated.resize(previous.size() + shift, 0);
        }
        for (size_t i = 0; i < previous.size(); ++i) {
            updated[i + shift] ^= gf_mul(scale, previous[i]);
        }
        if (2 * errors <= step) {
            previous = lambda;
            errors = step + 1 - errors;
            previous_discrepancy = discrepancy;
            shift = 1;
        } else {
            ++shift;
        }
        lambda = updated;
    }
    if (2 * errors > parity) {
        return false;
    }

    // error evaluator omega = S(x) lambda(x) mod x^parity
    std::vector<uint8_t> omega(parity, 0);
    for (size_t i = 0; i < parity; ++i) {
        for (size_t k = 0; k <= i && k < lambda.size(); ++k) {
            omega[i] ^= gf_mul(lambda[k], syndromes[i - k]);
        }
    }
    // formal derivative of lambda: odd terms only in characteristic 2
    std::vector<uint8_t> derivative(lambda.size() > 1 ? lambda.size() - 1 : 1, 0);
    for (size_t i = 1; i < lambda.size(); i += 2) {
        derivative[i - 1] = lambda[i];
    }

    // Chien search and Forney: the byte at index k has degree e = n - 1 - k, locator X = alpha^e
    size_t found = 0;
    for (int k = 0; k < n; ++k) {
        const int e = n - 1 - k;
        const uint8_t x_inverse = gf_pow_alpha(-e);
        if (eval_low_first(lambda, x_inverse) != 0) {
            continue;
        }
        const uint8_t denominator = eval_low_first(derivative, x_inverse);
        if (denominator == 0) {
            return false;
        }
        c[k] ^= gf_mul(gf_pow_alpha(e), gf_div(eval_low_first(omega, x_inverse), denominator));
        ++found;
    }
    return found == errors;
}

// ---------------------------------------------------------------- checks

static int failures = 0;

static void check(bool condition, const char* what) {
    std::printf("%s: %s\n", condition ? "ok  " : "FAIL", what);
    if (!condition) {
        ++failures;
    }
}

// multiplication in GF(256) the slow way: carry-less product reduced by 0x11D
static uint8_t reference_mul(uint8_t a, uint8_t b) {
    unsigned product = 0;
    for (int i = 0; i < 8; ++i) {
        if (b & (1u << i)) {
            product ^= static_cast<unsigned>(a) << i;
        }
    }
    for (int i = 15; i >= 8; --i) {
        if (product & (1u << i)) {
            product ^= 0x11Du << (i - 8);
        }
    }
    return static_cast<uint8_t>(product);
}

int main() {
    constexpr auto exp = EXP;

    // alpha^0 .. alpha^15, as in every QR code reference
    const uint8_t expected_head[16] = {1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38};
    bool head_ok = true;
    for (int i = 0; i < 16; ++i) {
        head_ok = head_ok && exp[i] == expected_head[i];
    }
    check(head_ok, "exp table starts 01 02 04 ... 80 1d 3a");

    // alpha generates the 255 nonzero elements
    std::array<bool, 256> seen{};
    for (uint8_t value : exp) {
        seen[value] = true;
    }
    bool all_nonzero = !seen[0];
    for (int v = 1; v < 256; ++v) {
        all_nonzero = all_nonzero && seen[v];
    }
    check(all_nonzero, "alpha generates the 255 nonzero elements");

    // the exp table agrees with carry-less multiplication: alpha^(i+1) = alpha^i * 2
    bool mul_ok = true;
    for (int i = 0; i < 254; ++i) {
        mul_ok = mul_ok && exp[i + 1] == reference_mul(exp[i], 2);
    }
    check(mul_ok, "exp table matches carry-less multiplication modulo 0x11D");

    // generator polynomial for 10 error correction codewords (QR version 1-M)
    constexpr auto g = generator_bytes<10>();
    const uint8_t expected_g[11] = {1, 216, 194, 159, 111, 199, 94, 95, 113, 157, 193};
    bool g_ok = true;
    for (int i = 0; i < 11; ++i) {
        g_ok = g_ok && g[i] == expected_g[i];
    }
    check(g_ok, "generator polynomial (x - a^0)...(x - a^9) computed at compile time");

    // "HELLO WORLD", QR code version 1-M: 16 data codewords, 10 error correction codewords
    const std::vector<uint8_t> data = {32, 91, 11, 120, 209, 114, 220, 77, 67, 64, 236, 17, 236, 17, 236, 17};
    const std::vector<uint8_t> ec = {196, 35, 39, 119, 235, 215, 231, 226, 93, 23};
    std::vector<uint8_t> codeword = rs_encode<10>(data);
    check(std::vector<uint8_t>(codeword.begin() + 16, codeword.end()) == ec,
          "encoding HELLO WORLD gives the error correction codewords of the QR standard");

    // 5 corrupted bytes (the most 10 parity bytes can correct)
    std::vector<uint8_t> received = codeword;
    const size_t positions[5] = {0, 7, 13, 20, 25};
    for (size_t i = 0; i < 5; ++i) {
        received[positions[i]] ^= static_cast<uint8_t>(17 * (i + 3));
    }
    const bool decoded = rs_decode<10>(&received);
    check(decoded && received == codeword, "5 corrupted bytes are corrected");

    // a clean codeword decodes unchanged
    std::vector<uint8_t> clean = codeword;
    check(rs_decode<10>(&clean) && clean == codeword, "a clean codeword decodes unchanged");

    // random messages with 0 to 5 corrupted bytes, at random positions
    std::mt19937 rng(2026);
    std::uniform_int_distribution<int> byte(0, 255), length(1, 200), count(0, 5);
    int corrected = 0;
    const int trials = 10000;
    for (int trial = 0; trial < trials; ++trial) {
        std::vector<uint8_t> message(length(rng));
        for (uint8_t& b : message) {
            b = static_cast<uint8_t>(byte(rng));
        }
        const std::vector<uint8_t> sent = rs_encode<10>(message);
        std::vector<uint8_t> noisy = sent;
        std::uniform_int_distribution<size_t> position(0, sent.size() - 1);
        for (int e = count(rng); e > 0; --e) {
            noisy[position(rng)] ^= static_cast<uint8_t>(byte(rng) | 1);
        }
        corrected += rs_decode<10>(&noisy) && noisy == sent ? 1 : 0;
    }
    check(corrected == trials, "10000 random codewords with up to 5 corrupted bytes are corrected");

    std::printf("%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
