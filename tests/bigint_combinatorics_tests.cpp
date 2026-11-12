// factorials, binomial coefficients and the other combinatorial numbers over bigint: exact beyond 64 bits.
// Expected values come from Python's integers (math.factorial, math.comb, fractions.Fraction).
#include <gtest/gtest.h>
#include <cstdint>
#include <string>
#include <type_traits>

#include "aerobus.h"

namespace {
using aerobus::bigint;
using aerobus::i64;

template<aerobus::string_literal S>
using N = bigint::from_string_t<S>;
}  // namespace

TEST(bigint_combinatorics, factorial) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // the type form is exact, whatever the size of the result
    static_assert(std::is_same_v<factorial_t<bigint, 0>, bigint::one>);
    static_assert(std::is_same_v<factorial_t<bigint, 1>, bigint::one>);
    static_assert(std::is_same_v<factorial_t<bigint, 20>, N<"2432902008176640000">>);
    // 21! is the first one that does not fit in an int64_t
    static_assert(std::is_same_v<factorial_t<bigint, 21>, N<"51090942171709440000">>);
    static_assert(std::is_same_v<factorial_t<bigint, 25>, N<"15511210043330985984000000">>);
    static_assert(std::is_same_v<factorial_t<bigint, 30>, N<"265252859812191058636308480000000">>);
    static_assert(std::is_same_v<factorial_t<bigint, 50>,
        N<"30414093201713378043612608166064768844377641568960512000000000000">>);

    EXPECT_EQ((factorial_t<bigint, 21>::to_string()), "51090942171709440000");
    EXPECT_EQ((factorial_t<bigint, 100>::to_string()),
        "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518"
        "286253697920827223758251185210916864000000000000000000000000");
}

TEST(bigint_combinatorics, factorial_agrees_with_native_integers) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // where the result fits in an int64_t, bigint and i64 compute the same thing
    static_assert(factorial_v<bigint, 20, int64_t> == factorial_v<i64, 20>);
    static_assert(combination_v<bigint, 30, 60, int64_t> == combination_v<i64, 30, 60>);
    static_assert(bell_v<bigint, 25, int64_t> == bell_v<i64, 25>);
    static_assert(stirling_2_v<bigint, 9, 4, int64_t> == stirling_2_v<i64, 9, 4>);
    static_assert(stirling_1_signed_v<bigint, 9, 4, int64_t> == stirling_1_signed_v<i64, 9, 4>);
    EXPECT_EQ((factorial_v<bigint, 20, int64_t>), 2432902008176640000LL);
}

TEST(bigint_combinatorics, factorial_as_value) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // bigint has no native value type: the result type is named, and conversion is correctly rounded
    static_assert(factorial_v<bigint, 20, int64_t> == 2432902008176640000LL);
    static_assert(factorial_v<bigint, 30, double> == 265252859812191058636308480000000.0);
    static_assert(factorial_v<bigint, 0, int> == 1);
    EXPECT_DOUBLE_EQ((factorial_v<bigint, 30, double>), 2.6525285981219107e32);
    EXPECT_FLOAT_EQ((factorial_v<bigint, 30, float>), 2.6525286e32f);
}

TEST(bigint_combinatorics, combination) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    static_assert(std::is_same_v<combination_t<bigint, 0, 10>, bigint::one>);
    static_assert(std::is_same_v<combination_t<bigint, 5, 10>, N<"252">>);
    static_assert(std::is_same_v<combination_t<bigint, 10, 40>, N<"847660528">>);
    static_assert(std::is_same_v<combination_t<bigint, 30, 60>, N<"118264581564861424">>);
    static_assert(std::is_same_v<combination_t<bigint, 50, 100>, N<"100891344545564193334812497256">>);
    // symmetric, and zero when there is nothing to choose from
    static_assert(std::is_same_v<combination_t<bigint, 50, 100>, combination_t<bigint, 50, 100>>);
    static_assert(std::is_same_v<combination_t<bigint, 90, 100>, combination_t<bigint, 10, 100>>);
    static_assert(std::is_same_v<combination_t<bigint, 5, 4>, bigint::zero>);
    static_assert(std::is_same_v<combination_t<bigint, 1, 0>, bigint::zero>);

    EXPECT_EQ((combination_t<bigint, 50, 100>::to_string()), "100891344545564193334812497256");
    static_assert(combination_v<bigint, 50, 100, double> == 100891344545564193334812497256.0);
    static_assert(combination_v<bigint, 5, 10, int> == 252);
}

TEST(bigint_combinatorics, pascal_rule) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // C(k, n) = C(k - 1, n - 1) + C(k, n - 1), on values that do not fit in an int64_t
    static_assert(std::is_same_v<
        combination_t<bigint, 50, 100>,
        bigint::add_t<combination_t<bigint, 49, 99>, combination_t<bigint, 50, 99>>>);
    // the sum of a row is a power of two
    static_assert(std::is_same_v<
        bigint::add_t<combination_t<bigint, 0, 2>,
            bigint::add_t<combination_t<bigint, 1, 2>, combination_t<bigint, 2, 2>>>,
        N<"4">>);
}

TEST(bigint_combinatorics, factorial_identities) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // n! = n * (n - 1)!, beyond 64 bits
    static_assert(std::is_same_v<
        factorial_t<bigint, 30>,
        bigint::mul_t<bigint::inject_constant_t<30>, factorial_t<bigint, 29>>>);
    // C(k, n) = n! / (k! (n - k)!), as an exact division
    static_assert(std::is_same_v<
        combination_t<bigint, 10, 40>,
        bigint::div_t<factorial_t<bigint, 40>,
            bigint::mul_t<factorial_t<bigint, 10>, factorial_t<bigint, 30>>>>);
    // and the division leaves no remainder
    static_assert(std::is_same_v<
        bigint::mod_t<factorial_t<bigint, 40>,
            bigint::mul_t<factorial_t<bigint, 10>, factorial_t<bigint, 30>>>,
        bigint::zero>);
}

TEST(bigint_combinatorics, inject_constant_from_size_t) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // the combinatorial templates index with size_t: an unsigned constant injects without a sign
    static_assert(std::is_same_v<bigint::inject_constant_t<static_cast<size_t>(0)>, bigint::zero>);
    static_assert(std::is_same_v<bigint::inject_constant_t<static_cast<size_t>(1)>, bigint::one>);
    static_assert(std::is_same_v<bigint::inject_constant_t<SIZE_MAX>,
        bigint::val<false, 0xffffffffU, 0xffffffffU>>);
    static_assert(std::is_same_v<bigint::inject_constant_t<static_cast<size_t>(0x100000000ULL)>,
        bigint::val<false, 0x1U, 0x0U>>);
    // signed constants keep their sign
    static_assert(std::is_same_v<bigint::inject_constant_t<-1>, bigint::val<true, 0x1U>>);
}

TEST(bigint_combinatorics, bell_numbers) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    static_assert(std::is_same_v<bell_t<bigint, 0>, bigint::one>);
    static_assert(std::is_same_v<bell_t<bigint, 1>, bigint::one>);
    static_assert(std::is_same_v<bell_t<bigint, 10>, N<"115975">>);
    static_assert(std::is_same_v<bell_t<bigint, 25>, N<"4638590332229999353">>);
    static_assert(std::is_same_v<bell_t<bigint, 30>, N<"846749014511809332450147">>);
    EXPECT_EQ((bell_t<bigint, 30>::to_string()), "846749014511809332450147");
    static_assert(bell_v<bigint, 10, int64_t> == 115975);
}

TEST(bigint_combinatorics, stirling_numbers) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    static_assert(std::is_same_v<stirling_2_t<bigint, 9, 4>, N<"7770">>);
    static_assert(std::is_same_v<stirling_2_t<bigint, 25, 5>, N<"2436684974110751">>);
    static_assert(std::is_same_v<stirling_2_t<bigint, 30, 7>, N<"4168916722553086402080">>);
    static_assert(std::is_same_v<stirling_1_signed_t<bigint, 5, 3>, N<"35">>);
    static_assert(std::is_same_v<stirling_1_signed_t<bigint, 25, 5>, N<"2677503356427960382362624">>);
    static_assert(std::is_same_v<stirling_1_signed_t<bigint, 30, 7>, N<"-11139316913434780466101123891200">>);
    static_assert(std::is_same_v<stirling_1_unsigned_t<bigint, 30, 7>, N<"11139316913434780466101123891200">>);
    EXPECT_EQ((stirling_1_signed_t<bigint, 30, 7>::to_string()), "-11139316913434780466101123891200");
    static_assert(stirling_2_v<bigint, 9, 4, int64_t> == 7770);
    static_assert(stirling_1_unsigned_v<bigint, 9, 4, int64_t> == 67284);
}

TEST(bigint_combinatorics, alternate_and_pow) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    static_assert(std::is_same_v<alternate_t<bigint, 0>, bigint::one>);
    static_assert(std::is_same_v<alternate_t<bigint, 1>, bigint::val<true, 0x1U>>);
    static_assert(std::is_same_v<alternate_t<bigint, 2>, bigint::one>);
    static_assert(alternate_v<bigint, 3, int> == -1);
    static_assert(std::is_same_v<pow_t<bigint::inject_constant_t<2>, 100>,
        N<"1267650600228229401496703205376">>);
    static_assert(std::is_same_v<pow_t<bigint::inject_constant_t<3>, 200>,
        N<"265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699044001">>);
    static_assert(pow_v<bigint::inject_constant_t<2>, 100, double> == 1267650600228229401496703205376.0);
}

TEST(bigint_combinatorics, bernoulli_numbers) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // B(n) lives in FractionField<bigint>, where numerator and denominator stay exact
    using B10 = bernoulli_t<bigint, 10>;
    static_assert(std::is_same_v<B10::x, N<"5">>);
    static_assert(std::is_same_v<B10::y, N<"66">>);
    using B20 = bernoulli_t<bigint, 20>;
    static_assert(std::is_same_v<B20::x, N<"-174611">>);
    static_assert(std::is_same_v<B20::y, N<"330">>);
    using B30 = bernoulli_t<bigint, 30>;
    static_assert(std::is_same_v<B30::x, N<"8615841276005">>);
    static_assert(std::is_same_v<B30::y, N<"14322">>);
    EXPECT_DOUBLE_EQ((bernoulli_v<double, bigint, 30>), 601580873.90064239);
}

TEST(bigint_combinatorics, taylor_coefficients_over_bigint) {
    using namespace aerobus;  // NOLINT(build/namespaces)
    // what the factorials are there for: a Taylor series whose coefficients do not fit in native integers
    using S = aerobus::sin<bigint, 25>;
    EXPECT_NEAR(S::eval(0.5), 0.479425538604203, 1e-15);
    using E = aerobus::exp<bigint, 25>;
    EXPECT_NEAR(E::eval(1.0), 2.718281828459045, 1e-14);
    // the coefficient of x^21 in exp is 1/21!, exactly
    using C21 = E::coeff_at_t<21>;
    static_assert(std::is_same_v<C21::x, bigint::one>);
    static_assert(std::is_same_v<C21::y, N<"51090942171709440000">>);
}
