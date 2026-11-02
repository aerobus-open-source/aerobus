#include <gtest/gtest.h>
#include <concepts>
#include <string>
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

namespace {
// every member a built-in ring provides (test-only: user rings need only the concepts)
template<typename R>
concept FullRingInterface = IsEuclideanDomain<R> && requires {
    typename R::template add_t<typename R::one, typename R::one>;
    typename R::template gt_t<typename R::one, typename R::zero>;
    typename R::template lt_t<typename R::one, typename R::zero>;
    typename R::template inject_constant_t<2>;
    typename R::one::enclosing_type;
    typename R::one::is_zero_t;
    requires std::is_same_v<typename R::one::enclosing_type, R>;
    { R::template eq_v<typename R::one, typename R::one> } -> std::convertible_to<bool>;
    { R::template gt_v<typename R::one, typename R::zero> } -> std::convertible_to<bool>;
    { R::template lt_v<typename R::one, typename R::zero> } -> std::convertible_to<bool>;
    { R::is_field } -> std::convertible_to<bool>;
    { R::is_integral_domain } -> std::convertible_to<bool>;
    { R::one::is_zero_v } -> std::convertible_to<bool>;
    { R::one::to_string() } -> std::convertible_to<std::string>;
};
}  // namespace

TEST(concepts, built_in_rings_share_one_interface) {
    EXPECT_TRUE(FullRingInterface<i16>);
    EXPECT_TRUE(FullRingInterface<i32>);
    EXPECT_TRUE(FullRingInterface<i64>);
    EXPECT_TRUE(FullRingInterface<zpz<7>>);
    EXPECT_TRUE(FullRingInterface<q32>);
    EXPECT_TRUE(FullRingInterface<q64>);
    EXPECT_TRUE(FullRingInterface<pi64>);
    EXPECT_TRUE(FullRingInterface<pq64>);
    EXPECT_TRUE((FullRingInterface<Quotient<i64, i64::val<17>>>));
    EXPECT_TRUE(FullRingInterface<FractionField<pq64>>);
}

TEST(concepts, zpz_injects_integers) {
    EXPECT_TRUE((std::is_same_v<zpz<7>::inject_ring_t<i32::val<9>>, zpz<7>::val<2>>));
}
