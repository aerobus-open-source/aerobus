#include <gtest/gtest.h>
#include <type_traits>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(concepts, ring) {
    EXPECT_TRUE((aerobus::IsRing<i32>));
    EXPECT_TRUE((aerobus::IsRing<i64>));
    EXPECT_TRUE((aerobus::IsRing<zpz<3>>));
    EXPECT_TRUE((aerobus::IsRing<aerobus::q32>));
    EXPECT_TRUE((aerobus::IsField<aerobus::q32>));
    EXPECT_TRUE((aerobus::IsRing<aerobus::polynomial<i32>>));
    EXPECT_TRUE((aerobus::IsField<aerobus::FractionField<aerobus::polynomial<i32>>>));
}

TEST(concepts, is_field_checks_the_flag) {
    EXPECT_FALSE((IsField<i32>));
    EXPECT_FALSE((IsField<i64>));
    EXPECT_FALSE((IsField<zpz<4>>));
    EXPECT_TRUE((IsField<zpz<5>>));
    EXPECT_FALSE((IsField<polynomial<i32>>));
    EXPECT_FALSE((IsField<polynomial<q32>>));
    EXPECT_TRUE((IsField<q64>));
}

namespace {
// has every member IsEuclideanDomain asks for; the flags are template parameters
template<bool euclidean, bool positive_one>
struct flagged_ring {
    template<int x>
    struct val {
        using enclosing_type = flagged_ring;
    };
    using zero = val<0>;
    using one = val<1>;
    template<typename, typename> using add_t = one;
    template<typename, typename> using sub_t = one;
    template<typename, typename> using mul_t = one;
    template<typename, typename> using div_t = one;
    template<typename, typename> using mod_t = zero;
    template<typename, typename> using gcd_t = one;
    template<typename, typename> using eq_t = std::true_type;
    template<typename> using pos_t = std::bool_constant<positive_one>;
    template<typename> static constexpr bool pos_v = positive_one;
    static constexpr bool is_euclidean_domain = euclidean;
    static constexpr bool is_field = false;
};
}  // namespace

TEST(concepts, is_euclidean_domain_checks_the_flags) {
    EXPECT_TRUE((IsRing<flagged_ring<false, true>>));
    EXPECT_TRUE((IsEuclideanDomain<flagged_ring<true, true>>));
    EXPECT_FALSE((IsEuclideanDomain<flagged_ring<false, true>>));
    EXPECT_FALSE((IsEuclideanDomain<flagged_ring<true, false>>));
}
