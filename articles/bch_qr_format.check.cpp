// Appended to the code blocks of bch_qr_format.md by tools/check_article.py:
// the numbers and the table quoted in the article are checked by the compiler.

template<int32_t... coeffs>
using P2 = aerobus::make_int_polynomial_t<F2, coeffs...>;

static_assert(std::is_same_v<factors, aerobus::type_list<
    aerobus::factor<P2<1, 1>, 1>, aerobus::factor<P2<1, 1, 1>, 1>, aerobus::factor<P2<1, 0, 0, 1, 1>, 1>,
    aerobus::factor<P2<1, 1, 0, 0, 1>, 1>, aerobus::factor<P2<1, 1, 1, 1, 1>, 1>>>,
    "the five factors of x^15 - 1, in the order of the table");

template<typename f, size_t k>
constexpr bool root = vanishes_at<f, k>(std::make_index_sequence<f::degree + 1>{});

static_assert(root<P2<1, 1>, 0> && root<P2<1, 1>, 15>);
static_assert(root<P2<1, 1, 1>, 5> && root<P2<1, 1, 1>, 10>);
static_assert(root<P2<1, 0, 0, 1, 1>, 1> && root<P2<1, 0, 0, 1, 1>, 2> && root<P2<1, 0, 0, 1, 1>, 4> &&
              root<P2<1, 0, 0, 1, 1>, 8>);
static_assert(root<P2<1, 1, 0, 0, 1>, 7> && root<P2<1, 1, 0, 0, 1>, 14> && root<P2<1, 1, 0, 0, 1>, 13> &&
              root<P2<1, 1, 0, 0, 1>, 11>);
static_assert(root<P2<1, 1, 1, 1, 1>, 3> && root<P2<1, 1, 1, 1, 1>, 6> && root<P2<1, 1, 1, 1, 1>, 12> &&
              root<P2<1, 1, 1, 1, 1>, 9>);

// "changing 6 to 4 keeps only two factors"
template<typename... Fs>
constexpr size_t kept(aerobus::type_list<Fs...>, auto roots) {
    return ((has_root_among<typename Fs::type>(roots) ? size_t{1} : size_t{0}) + ...);
}
static_assert(kept(factors{}, std::make_index_sequence<6>{}) == 3);
static_assert(kept(factors{}, std::make_index_sequence<4>{}) == 2);

static_assert(G == 0x537);

// the table of ISO/IEC 18004, indexed by (level bits << 3) | mask: M, L, H, Q
constexpr std::array<uint16_t, 32> STANDARD = {
    0x5412, 0x5125, 0x5e7c, 0x5b4b, 0x45f9, 0x40ce, 0x4f97, 0x4aa0,
    0x77c4, 0x72f3, 0x7daa, 0x789d, 0x662f, 0x6318, 0x6c41, 0x6976,
    0x1689, 0x13be, 0x1ce7, 0x19d0, 0x0762, 0x0255, 0x0d0c, 0x083b,
    0x355f, 0x3068, 0x3f31, 0x3a06, 0x24b4, 0x2183, 0x2eda, 0x2bed,
};

constexpr bool words_match() {
    for (uint32_t data = 0; data < 32; ++data) {
        if (encode_format(data) != STANDARD[data]) {
            return false;
        }
    }
    return true;
}
static_assert(words_match(), "the 32 format words are those of the standard");
static_assert(encode_format(0) == 0x5412 && encode_format(8) == 0x77c4, "M and L with mask 0, as quoted");

constexpr int minimum_distance() {
    int best = 15;
    for (uint32_t i = 0; i < 32; ++i) {
        for (uint32_t j = i + 1; j < 32; ++j) {
            int d = 0;
            for (uint32_t v = encode_format(i) ^ encode_format(j); v != 0; v >>= 1) {
                d += static_cast<int>(v & 1U);
            }
            best = d < best ? d : best;
        }
    }
    return best;
}
static_assert(minimum_distance() == 7, "minimum distance 7, as quoted");
static_assert(1 + 15 + 105 + 455 == 576, "patterns of at most 3 flipped bits among 15");
