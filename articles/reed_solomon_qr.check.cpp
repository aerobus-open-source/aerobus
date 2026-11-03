// Appended to the code blocks of reed_solomon_qr.md by tools/check_article.py:
// the numbers quoted in the article are checked by the compiler.

template<typename G, size_t... i>
constexpr std::array<uint8_t, sizeof...(i)> generator_bytes(std::index_sequence<i...>) {
    return {pack<typename G::template coeff_at_t<G::degree - i>>(std::make_index_sequence<8>{})...};
}

constexpr auto G10 = generator_bytes<rs_generator<10>::type>(std::make_index_sequence<11>{});
static_assert(G10 == std::array<uint8_t, 11>{1, 216, 194, 159, 111, 199, 94, 95, 113, 157, 193},
              "generator for 10 parity bytes, as quoted in the article");
static_assert(EXP[0] == 0x01 && EXP[1] == 0x02 && EXP[7] == 0x80 && EXP[9] == 0x3a,
              "exp table starts 01 02 04 ... 80 1d 3a");
