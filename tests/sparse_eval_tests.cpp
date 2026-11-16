// sparse_eval: P(x) = x^k Q(x^n) evaluated as x^k Q(x^n) (issue #23)
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

namespace {
using PZ = polynomial<i64>;
template<int64_t... cs>
using P = make_int_polynomial_t<i64, cs...>;

// x^6 + x^3 + 1 = Q(x^3) with Q = x^2 + x + 1
using P631 = P<1, 0, 0, 1, 0, 0, 1>;
// x^7 + x^3 = x^3 Q(x^4) with Q = x + 1
using P73 = P<1, 0, 0, 0, 1, 0, 0, 0>;
// 3x^2 + 2x + 1: dense
using Pdense = P<3, 2, 1>;
// 5x^4: a monomial
using Pmono = P<5, 0, 0, 0, 0>;
using Pconst = P<5>;
using Pzero = PZ::zero;

using Sin = aerobus::sin<i64, 9>;
using Cos = aerobus::cos<i64, 8>;
using Atan = aerobus::atan<i64, 9>;
using Cosh = aerobus::cosh<i64, 8>;
}  // namespace

TEST(sparse_eval, valuation_and_stride) {
    static_assert(P631::valuation == 0 && P631::stride == 3);
    static_assert(P73::valuation == 3 && P73::stride == 4);
    static_assert(Pdense::valuation == 0 && Pdense::stride == 1);
    static_assert(Pmono::valuation == 4 && Pmono::stride == 1);
    static_assert(Pconst::valuation == 0 && Pconst::stride == 1);
    static_assert(Pzero::valuation == 0 && Pzero::stride == 1);
    static_assert(Sin::valuation == 1 && Sin::stride == 2);
    static_assert(Cos::valuation == 0 && Cos::stride == 2);
    static_assert(Atan::valuation == 1 && Atan::stride == 2);
    static_assert(Cosh::valuation == 0 && Cosh::stride == 2);
    // x^2 + x^0: k = 0, gcd(2) = 2
    static_assert(P<1, 0, 1>::stride == 2);
    // x^5 + x^2: k = 2, gcd(3) = 3
    static_assert(P<1, 0, 0, 1, 0, 0>::valuation == 2 && P<1, 0, 0, 1, 0, 0>::stride == 3);
    // x^4 + x^2 + x: k = 1, gcd(3, 1) = 1
    static_assert(P<1, 0, 1, 1, 0>::valuation == 1 && P<1, 0, 1, 1, 0>::stride == 1);
}

TEST(sparse_eval, stride_reduced) {
    static_assert(std::is_same_v<P631::stride_reduced_t, P<1, 1, 1>>);
    static_assert(std::is_same_v<P73::stride_reduced_t, P<1, 1>>);
    static_assert(std::is_same_v<Pdense::stride_reduced_t, Pdense>);
    static_assert(std::is_same_v<Pmono::stride_reduced_t, P<5>>);
    static_assert(std::is_same_v<Pconst::stride_reduced_t, Pconst>);
    static_assert(std::is_same_v<Pzero::stride_reduced_t, Pzero>);
    using Q = polynomial<q64>;
    static_assert(std::is_same_v<Sin::stride_reduced_t,
        Q::val<Sin::coeff_at_t<9>, Sin::coeff_at_t<7>, Sin::coeff_at_t<5>, Sin::coeff_at_t<3>, Sin::coeff_at_t<1>>>);
    static_assert(std::is_same_v<Cos::stride_reduced_t,
        Q::val<Cos::coeff_at_t<8>, Cos::coeff_at_t<6>, Cos::coeff_at_t<4>, Cos::coeff_at_t<2>, Cos::coeff_at_t<0>>>);
}

TEST(sparse_eval, constexpr_and_exact_on_integers) {
    static_assert(P73::sparse_eval<int64_t>(2) == 136);
    static_assert(P631::sparse_eval<int64_t>(2) == 73);
    static_assert(Pmono::sparse_eval<int64_t>(3) == 405);
    static_assert(Pconst::sparse_eval<int64_t>(3) == 5);
    static_assert(Pzero::sparse_eval<int64_t>(3) == 0);
    for (int64_t x = -20; x <= 20; ++x) {
        EXPECT_EQ(P73::sparse_eval(x), P73::eval(x));
        EXPECT_EQ(P631::sparse_eval(x), P631::eval(x));
        EXPECT_EQ(Pdense::sparse_eval(x), Pdense::eval(x));
        EXPECT_EQ(Pmono::sparse_eval(x), Pmono::eval(x));
        EXPECT_EQ((P<1, 0, 0, 1, 0, 0>::sparse_eval(x)), (P<1, 0, 0, 1, 0, 0>::eval(x)));
    }
    EXPECT_EQ((P<1, 0, 0>::sparse_eval<int16_t>(100)), 10000);
}

TEST(sparse_eval, dense_is_eval) {
    for (int i = -100; i <= 100; ++i) {
        const double x = i / 37.0;
        EXPECT_EQ(Pdense::sparse_eval(x), Pdense::eval(x));
    }
}

namespace {
// error of approx in units of the last place of ref, rounded to T
template<typename T>
double ulps(T approx, long double ref) {
    const T r = static_cast<T>(ref);
    const T ulp = std::nextafter(std::fabs(r), std::numeric_limits<T>::infinity()) - std::fabs(r);
    return static_cast<double>(std::fabs(static_cast<long double>(approx) - ref) / ulp);
}

// about as accurate as Horner: a zero coefficient costs Horner a rounded product x * accum, and
// sparse_eval rounds x^n once instead, so the maximum error moves by a fraction of an ulp either way
template<typename Poly, typename T>
void expect_as_accurate(const char* name) {
    double max_horner = 0, max_sparse = 0;
    constexpr int N = 20000;
    for (int i = -N; i <= N; ++i) {
        const T x = static_cast<T>(i) / static_cast<T>(N);
        const long double ref = Poly::template eval<long double>(static_cast<long double>(x));
        if (ref == 0) {
            continue;
        }
        max_horner = std::max(max_horner, ulps<T>(Poly::eval(x), ref));
        max_sparse = std::max(max_sparse, ulps<T>(Poly::sparse_eval(x), ref));
    }
    EXPECT_LE(max_sparse, max_horner + 0.25) << name;
    EXPECT_LT(max_sparse, 2.0) << name;
}
}  // namespace

TEST(sparse_eval, as_accurate_as_horner_float) {
    expect_as_accurate<Sin, float>("sin");
    expect_as_accurate<Cos, float>("cos");
    expect_as_accurate<Atan, float>("atan");
    expect_as_accurate<Cosh, float>("cosh");
}

TEST(sparse_eval, as_accurate_as_horner_double) {
    if constexpr (std::numeric_limits<long double>::digits < 64) {
        GTEST_SKIP() << "needs a long double with at least 64 significant bits";
    } else {
        expect_as_accurate<Sin, double>("sin");
        expect_as_accurate<Cos, double>("cos");
        expect_as_accurate<Atan, double>("atan");
        expect_as_accurate<Cosh, double>("cosh");
    }
}
