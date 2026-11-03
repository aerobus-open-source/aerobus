#include <gtest/gtest.h>
#include <type_traits>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(embedding, simple) {
    {
        using Small = Quotient<i32, i32::val<4>>;
        using embedded = typename Embed<Small, i32>::type<typename Small::template inject_constant_t<3>>;
        EXPECT_TRUE((std::is_same_v<embedded, typename i32::template inject_constant_t<3>>));
    }
    {
        using Small = zpz<4>;
        using embedded = typename Embed<Small, i32>::type<typename Small::template inject_constant_t<3>>;
        EXPECT_EQ(embedded::v, 3);
    }
    {
        using embedded = typename Embed<i32, i64>::type<i32::val<12>>;
        EXPECT_TRUE((std::is_same_v<embedded, typename i64::template inject_constant_t<12>>));
    }
    {
        using embedded = typename Embed<i32, q32>::type<i32::val<12>>;
        EXPECT_TRUE((std::is_same_v<embedded, q32::inject_constant_t<12>>));
    }
    {
        using embedded = typename Embed<q32, q64>::type<make_q32_t<1, 2>>;
        EXPECT_TRUE((std::is_same_v<embedded, make_q64_t<1, 2>>));
    }
}

TEST(embedding, polynomial) {
    {  // 1+X
        using embedded = Embed<polynomial<i32>, polynomial<i64>>::type<polynomial<i32>::val<i32::one, i32::one>>;
        EXPECT_TRUE((std::is_same_v<embedded, polynomial<i64>::val<i64::one, i64::one>>));
    }
    {  // 1+2X
        using embedded = Embed<polynomial<i32>, polynomial<i64>>::type<polynomial<i32>::val<
                            i32::inject_constant_t<2>, i32::one>>;
        EXPECT_TRUE((std::is_same_v<embedded, polynomial<i64>::val<i64::inject_constant_t<2>, i64::one>>));
    }
    {  // 1+2X, from i32 to q32
        using embedded = embed_int_poly_in_fractions_t<make_int_polynomial_t<i32, 2, 1>>;
        EXPECT_TRUE((std::is_same_v<embedded, make_frac_polynomial_t<i32, 2, 1>>))
            << "actual : " << embedded::to_string();
    }
}

TEST(embedding, quotient_embeds_its_canonical_representative) {
    using Z17 = Quotient<i64, i64::val<17>>;
    using sum = Z17::add_t<Z17::inject_constant_t<16>, Z17::inject_constant_t<16>>;
    EXPECT_TRUE((std::is_same_v<Embed<Z17, i64>::type<sum>, i64::val<15>>));
}
