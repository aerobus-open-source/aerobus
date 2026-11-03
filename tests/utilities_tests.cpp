#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(utilities, type_at) {
    EXPECT_TRUE((std::is_same_v<internal::type_at_t<0, float, int, int64_t>, float>));
    EXPECT_TRUE((std::is_same_v<internal::type_at_t<1, float, int, int64_t>, int>));
    EXPECT_TRUE((std::is_same_v<internal::type_at_t<2, float, int, int64_t>, int64_t>));
}

TEST(utilities, factorial) {
    constexpr float x = factorial_v<i32, 3>;
    EXPECT_EQ(x, 6.0f);

    constexpr size_t y = factorial_v<i32, 0>;
    EXPECT_EQ(y, 1);
}

TEST(utilities, combination) {
    constexpr int x = combination_v<i32, 2, 4>;
    EXPECT_EQ(x, 6);

    constexpr int y = combination_v<i32, 0, 4>;
    EXPECT_EQ(y, 1);

    constexpr int z = combination_v<i32, 1, 4>;
    EXPECT_EQ(z, 4);

    constexpr int zz = combination_v<i32, 3, 4>;
    EXPECT_EQ(zz, 4);
}

TEST(utilities, bernoulli) {
    constexpr float b0 = bernoulli_v<float, i32, 0>;
    EXPECT_EQ(b0, 1.0f);

    constexpr float b1 = bernoulli_v<float, i32, 1>;
    EXPECT_EQ(b1, -0.5f);

    using B2 = bernoulli_t<i32, 2>;
    EXPECT_EQ(B2::x::v, 1);
    EXPECT_EQ(B2::y::v, 6);

    constexpr double b3 = bernoulli_v<double, i32, 3>;
    EXPECT_EQ(b3, 0.0);

    using B4 = bernoulli_t<i64, 4>;
    EXPECT_EQ(B4::x::v, -1);
    EXPECT_EQ(B4::y::v, 30);
}

TEST(utilities, bell) {
    constexpr int64_t B0 = bell_v<i64, 0>;
    EXPECT_EQ(B0, 1);
    constexpr int64_t B1 = bell_v<i64, 1>;
    EXPECT_EQ(B1, 1);
    constexpr int64_t B2 = bell_v<i64, 2>;
    EXPECT_EQ(B2, 2);
    constexpr int64_t B3 = bell_v<i64, 3>;
    EXPECT_EQ(B3, 5);
    constexpr int64_t B4 = bell_v<i64, 4>;
    EXPECT_EQ(B4, 15);
}

TEST(utilities, is_prime) {
    // 64-bit numbers, beyond trial division: Mersenne prime 2^61 - 1, the largest prime below 2^64,
    // a strong pseudoprime to every prime base up to 31, a Carmichael number and a product of two primes
    static_assert(is_prime_v<2305843009213693951ULL>);
    static_assert(is_prime_v<18446744073709551557ULL>);
    static_assert(!is_prime_v<3825123056546413051ULL>);
    static_assert(!is_prime_v<561>);
    static_assert(!is_prime_v<998244359987710471ULL>);
    EXPECT_TRUE(is_prime_v<2>);
    EXPECT_TRUE(is_prime_v<3>);
    EXPECT_TRUE(is_prime_v<7>);
    EXPECT_TRUE(is_prime_v<11>);
    EXPECT_TRUE(is_prime_v<13>);
    EXPECT_TRUE(is_prime_v<31>);
    EXPECT_TRUE(is_prime_v<7883>);
    EXPECT_TRUE(is_prime_v<7919>);
    EXPECT_TRUE(is_prime_v<7927>);
    EXPECT_TRUE(is_prime_v<1000423>);
    EXPECT_FALSE(is_prime_v<1>);
    EXPECT_FALSE(is_prime_v<4>);
    EXPECT_FALSE(is_prime_v<7884>);
    EXPECT_FALSE(is_prime_v<7928>);
}

TEST(utilities, is_prime_odd_composites) {
    // products of two odd primes and squares of primes
    EXPECT_FALSE(is_prime_v<9>);
    EXPECT_FALSE(is_prime_v<15>);
    EXPECT_FALSE(is_prime_v<25>);
    EXPECT_FALSE(is_prime_v<35>);
    EXPECT_FALSE(is_prime_v<49>);
    EXPECT_FALSE(is_prime_v<55>);
    EXPECT_FALSE(is_prime_v<77>);
    EXPECT_FALSE(is_prime_v<121>);
    EXPECT_FALSE(is_prime_v<143>);
    EXPECT_FALSE(is_prime_v<187>);
    EXPECT_FALSE(is_prime_v<62773913>);  // 7919 * 7927
    EXPECT_TRUE(is_prime_v<997>);
    EXPECT_TRUE(is_prime_v<65537>);
    EXPECT_TRUE(is_prime_v<2147483647>);
}

template<size_t... Is>
constexpr std::array<bool, sizeof...(Is)> is_prime_table(std::index_sequence<Is...>) {
    return {is_prime_v<Is>...};
}

TEST(utilities, is_prime_matches_trial_division_below_1000) {
    constexpr auto table = is_prime_table(std::make_index_sequence<1000>{});
    for (size_t n = 0; n < table.size(); ++n) {
        bool expected = n >= 2;
        for (size_t d = 2; d * d <= n; ++d) {
            if (n % d == 0) {
                expected = false;
                break;
            }
        }
        EXPECT_EQ(table[n], expected) << "n = " << n;
    }
}

TEST(utilities, alternate) {
    constexpr int a0 = alternate_v<i32, 0>;
    EXPECT_EQ(a0, 1);
    constexpr int a1 = alternate_v<i32, 1>;
    EXPECT_EQ(a1, -1);
    constexpr int a2 = alternate_v<i32, 2>;
    EXPECT_EQ(a2, 1);
}

TEST(utilities, pow) {
    constexpr int32_t a0 = pow_t<typename i32::inject_constant_t<2>, 3>::v;
    EXPECT_EQ(a0, 8);
    constexpr int32_t a1 = pow_t<typename i32::inject_constant_t<2>, 4>::v;
    EXPECT_EQ(a1, 16);
    constexpr int32_t a2 = pow_t<i32::zero, 0>::v;
    EXPECT_EQ(a2, 1);
    constexpr int64_t a3 = pow_v<typename i64::inject_constant_t<3>, 3>;
    EXPECT_EQ(a3, 27);
}

TEST(utilities, pow_scalar) {
    constexpr int32_t a0 = pow_scalar<int32_t, 2>(2);
    EXPECT_EQ(a0, 4);

    constexpr double x = pow_scalar<double, 3>(3.0);
    EXPECT_EQ(x, 27.0);

    constexpr double one = pow_scalar<size_t, 0>(12);
    EXPECT_EQ(one, 1);

    constexpr double zone = pow_scalar<size_t, 0>(0);
    EXPECT_EQ(zone, 1);
}

TEST(utilities, abs) {
    using a0 = abs_t<i32::inject_constant_t<-1>>;
    EXPECT_EQ(a0::v, 1);
    using a1 = abs_t<i64::inject_constant_t<1>>;
    EXPECT_EQ(a1::v, 1);
    using a2 = abs_t<pi64::val<i64::val<-1>, i64::val<1>>>;
    using expected = pi64::val<i64::val<1>, i64::val<-1>>;
    EXPECT_TRUE((std::is_same_v<a2, expected>));
}

TEST(utilities, stirling_first_kind) {
    constexpr int s00 = stirling_1_signed_v<i64, 0, 0>;
    EXPECT_EQ(s00, 1);
    constexpr int s01 = stirling_1_signed_v<i64, 0, 1>;
    EXPECT_EQ(s01, 0);
    constexpr int s11 = stirling_1_signed_v<i64, 1, 1>;
    EXPECT_EQ(s11, 1);
    constexpr int s53 = stirling_1_signed_v<i64, 5, 3>;
    EXPECT_EQ(s53, 35);
    constexpr int s94 = stirling_1_signed_v<i64, 9, 4>;
    EXPECT_EQ(s94, -67284);
    constexpr int u94 = stirling_1_unsigned_v<i64, 9, 4>;
    EXPECT_EQ(u94, 67284);
}

TEST(utilities, stirling_second_kind) {
    constexpr int s00 = stirling_2_v<i64, 0, 0>;
    EXPECT_EQ(s00, 1);
    constexpr int s01 = stirling_2_v<i64, 0, 1>;
    EXPECT_EQ(s01, 0);
    constexpr int s11 = stirling_2_v<i64, 1, 1>;
    EXPECT_EQ(s11, 1);
    constexpr int s21 = stirling_2_v<i64, 2, 1>;
    EXPECT_EQ(s21, 1);
    constexpr int s32 = stirling_2_v<i64, 3, 2>;
    EXPECT_EQ(s32, 3);
    constexpr int s53 = stirling_2_v<i64, 5, 3>;
    EXPECT_EQ(s53, 25);
    constexpr int s94 = stirling_2_v<i64, 9, 4>;
    EXPECT_EQ(s94, 7770);
}

TEST(continuous_fractions, basic_assertions) {
    // A001203
    constexpr double A_PI = PI_fraction::val;
    EXPECT_FALSE((::fabs(A_PI - M_PI) > 0));

    // A003417
    constexpr double A_E = E_fraction::val;
    EXPECT_FALSE((::fabs(A_E - M_E) > 0));

    constexpr double A_SQRT2 = SQRT2_fraction::val;
    EXPECT_FALSE((::fabs(A_SQRT2 - M_SQRT2) > 0));

    constexpr double A_SQRT3 = SQRT3_fraction::val;
    EXPECT_FALSE((::fabs(A_SQRT3 - 1.7320508075688772935) > 0));
}

TEST(utilities, counts_with_k_greater_than_n_are_zero) {
    EXPECT_EQ((combination_v<i32, 5, 4>), 0);
    EXPECT_EQ((combination_v<i64, 1, 0>), 0);
    EXPECT_EQ((stirling_2_v<i64, 1, 2>), 0);
    EXPECT_EQ((stirling_2_v<i64, 3, 7>), 0);
    EXPECT_EQ((stirling_1_signed_v<i64, 2, 5>), 0);
}

TEST(continuous_fractions, get_in_any_type) {
    using half_plus_one = ContinuedFraction<1, 2>;  // 1 + 1/2
    EXPECT_EQ((half_plus_one::get<float>()), 1.5f);
    EXPECT_EQ((half_plus_one::get<double>()), 1.5);
    EXPECT_EQ((ContinuedFraction<3>::get<float>()), 3.0f);
    EXPECT_EQ(SQRT2_fraction::get<float>(), static_cast<float>(SQRT2_fraction::val));
}

const void* bell_v_address_elsewhere() {
    return &bell_v<i64, 4>;
}

const void* stirling_2_v_address_elsewhere() {
    return &stirling_2_v<i64, 5, 3>;
}

TEST(utilities, pow_scalar_odd_and_negative) {
    EXPECT_EQ((pow_scalar<int, 3>(-2)), -8);
    EXPECT_EQ((pow_scalar<double, 10>(2.0)), 1024.0);
    EXPECT_EQ((pow_scalar<double, 1>(5.0)), 5.0);
}
