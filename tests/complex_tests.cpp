#include <gtest/gtest.h>
#include <complex>
#include <string>
#include <type_traits>
#include <utility>

#define AEROBUS_CONWAY_IMPORTS
#include "aerobus.h"

using namespace aerobus;  // NOLINT

namespace {
// the Gaussian integers Z[i]
using Zi = Complex<i32>;
template<int32_t a, int32_t b>
using z = Zi::val<i32::val<a>, i32::val<b>>;

// Q[i], a field
using Qi = Complex<q32>;
template<int64_t a, int64_t b, int64_t c, int64_t d>
using q = Qi::val<make_q32_t<a, b>, make_q32_t<c, d>>;

// F_p[i]
template<int32_t p>
using Ci = Complex<zpz<p>>;
template<int32_t p, int32_t a, int32_t b>
using f = typename Ci<p>::template val<typename zpz<p>::template val<a>, typename zpz<p>::template val<b>>;

template<typename V, typename T>
concept has_get = requires { V::template get<T>(); };

// every nonzero element u of F_p[i] satisfies u * (1 / u) = 1
template<int32_t p, size_t... k>
constexpr bool all_invertible(std::index_sequence<k...>) {
    using C = Ci<p>;
    return ((std::is_same_v<typename C::template mul_t<f<p, (k + 1) / p, (k + 1) % p>,
                                                      typename C::template div_t<typename C::one,
                                                                                 f<p, (k + 1) / p, (k + 1) % p>>>,
                            typename C::one>) && ...);
}
}  // namespace

TEST(complex, gaussian_integers) {
    static_assert(std::is_same_v<Zi::mul_t<Zi::i, Zi::i>, z<-1, 0>>);
    static_assert(std::is_same_v<Zi::add_t<z<1, 2>, z<3, -5>>, z<4, -3>>);
    static_assert(std::is_same_v<Zi::sub_t<z<1, 2>, z<3, -5>>, z<-2, 7>>);
    // (1 + 2i)(3 - 5i) = 3 - 5i + 6i - 10i² = 13 + i
    static_assert(std::is_same_v<Zi::mul_t<z<1, 2>, z<3, -5>>, z<13, 1>>);
    static_assert(std::is_same_v<Zi::conj_t<z<1, 2>>, z<1, -2>>);
    static_assert(std::is_same_v<Zi::norm_t<z<3, 4>>, i32::val<25>>);
    static_assert(Zi::eq_v<z<1, 2>, z<1, 2>> && !Zi::eq_v<z<1, 2>, z<2, 1>>);
    static_assert(Zi::val<i32::val<1>, i32::val<2>>::x::v == 1 && z<1, 2>::y::v == 2);
    static_assert(z<0, 0>::is_zero_t::value && z<0, 0>::is_zero_v && !Zi::i::is_zero_v);
    static_assert(std::is_same_v<Zi::zero, z<0, 0>> && std::is_same_v<Zi::one, z<1, 0>>);
    static_assert(std::is_same_v<Zi::inject_constant_t<3>, z<3, 0>>);
    static_assert(std::is_same_v<Zi::inject_ring_t<i32::val<-3>>, z<-3, 0>>);
    static_assert(std::is_same_v<Embed<i32, Zi>::type<i32::val<3>>, z<3, 0>>);
    static_assert(IsRing<Zi> && !IsEuclideanDomain<Zi>);
    static_assert(Zi::is_integral_domain && !Zi::is_field);
    static_assert(Complex<bigint>::is_integral_domain);
    // the generic operators
    static_assert(std::is_same_v<mul_t<Zi::i, Zi::i>, z<-1, 0>>);
}

TEST(complex, rationals_form_a_field) {
    static_assert(IsField<Qi>);
    static_assert(Complex<FractionField<bigint>>::is_field);
    // (1 + 2i) / (3 + 4i) = (1 + 2i)(3 - 4i) / 25 = (11 + 2i) / 25
    static_assert(std::is_same_v<Qi::div_t<q<1, 1, 2, 1>, q<3, 1, 4, 1>>, q<11, 25, 2, 25>>);
    static_assert(std::is_same_v<Qi::div_t<Qi::one, Qi::i>, q<0, 1, -1, 1>>);  // 1 / i = -i
    static_assert(std::is_same_v<Qi::mod_t<q<1, 1, 2, 1>, q<3, 1, 4, 1>>, Qi::zero>);
    static_assert(std::is_same_v<Qi::gcd_t<q<1, 1, 2, 1>, Qi::zero>, Qi::one>);
    static_assert(std::is_same_v<Qi::gcd_t<Qi::zero, Qi::zero>, Qi::zero>);
    static_assert(Qi::pos_v<Qi::i> && !Qi::pos_v<Qi::zero>);
    static_assert(std::is_same_v<divide_t<Qi::i, Qi::i>, Qi::one>);
    static_assert(std::is_same_v<pow_t<q<1, 1, 1, 1>, 4>, q<-4, 1, 0, 1>>);  // (1 + i)² = 2i
}

TEST(complex, finite_fields) {
    // x² + 1 is irreducible over F_p, and F_p[i] a field with p² elements, if and only if p = 3 mod 4
    static_assert(IsField<Ci<3>> && IsField<Ci<7>>);
    static_assert(all_invertible<3>(std::make_index_sequence<8>{}));
    static_assert(all_invertible<7>(std::make_index_sequence<48>{}));
    // in F_5, -1 = 2²: (2 + i)(2 - i) = 4 + 1 = 0
    static_assert(std::is_same_v<Ci<5>::mul_t<f<5, 2, 1>, f<5, 2, -1>>, Ci<5>::zero>);
    static_assert(!Ci<5>::is_field && !Ci<5>::is_integral_domain && IsRing<Ci<5>> && !IsEuclideanDomain<Ci<5>>);
    // in F_2, (1 + i)² = 2i = 0
    static_assert(std::is_same_v<Ci<2>::mul_t<f<2, 1, 1>, f<2, 1, 1>>, Ci<2>::zero>);
    static_assert(!Ci<2>::is_field && !Ci<2>::is_integral_domain);
    // the same rule over F_q: q = 27 = 3 mod 4, q = 9 = 1 mod 4
    static_assert(Complex<GaloisField<3, 3>>::is_field);
    static_assert(!Complex<GaloisField<3, 2>>::is_field);
    // -1 = i² is a square in Q[i]
    static_assert(!Complex<Qi>::is_integral_domain && !Complex<Qi>::is_field);
}

TEST(complex, polynomials) {
    using P = polynomial<Qi>;
    using X_minus_i = P::sub_t<P::X, P::val<Qi::i>>;
    using X_plus_i = P::add_t<P::X, P::val<Qi::i>>;
    using X2_plus_1 = P::val<Qi::one, Qi::zero, Qi::one>;
    static_assert(std::is_same_v<P::mul_t<X_minus_i, X_plus_i>, X2_plus_1>);
    static_assert(std::is_same_v<P::div_t<X2_plus_1, X_minus_i>, X_plus_i>);
    static_assert(std::is_same_v<P::mod_t<X2_plus_1, X_minus_i>, P::zero>);
}

TEST(complex, bigint_and_its_fractions) {
    using B = bigint;
    using ZB = Complex<B>;
    using QB = FractionField<B>;
    using CQ = Complex<QB>;
    using two_64 = B::val<false, 1U, 0U, 0U>;
    using two_128 = B::val<false, 1U, 0U, 0U, 0U, 0U>;
    using three = B::val<false, 3U>;
    using nine = B::val<false, 9U>;
    // a = 2^64 + 3i, a² = 2^128 - 9 + 3 · 2^65 i, |a|² = 2^128 + 9
    using a = ZB::val<two_64, three>;
    static_assert(std::is_same_v<ZB::mul_t<a, a>, ZB::val<B::sub_t<two_128, nine>, B::val<false, 6U, 0U, 0U>>>);
    static_assert(std::is_same_v<ZB::norm_t<a>, B::add_t<two_128, nine>>);
    static_assert(IsRing<ZB> && ZB::is_integral_domain && !ZB::is_field);
    // 1 / a = conj(a) / |a|², correctly rounded part by part
    using qa = CQ::val<QB::val<two_64, B::one>, QB::val<three, B::one>>;
    using inverse = CQ::div_t<CQ::one, qa>;
    static_assert(IsField<CQ>);
    static_assert(std::is_same_v<CQ::mul_t<inverse, qa>, CQ::one>);
    static_assert(std::is_same_v<inverse, CQ::val<QB::val<two_64, B::add_t<two_128, nine>>,
                                                  QB::val<B::val<true, 3U>, B::add_t<two_128, nine>>>>);
    static_assert(inverse::get<std::complex<double>>() ==
                  std::complex<double>(0x1p-64, -0x1.8p-127));
    // (X - a)(X - conj(a)) has real coefficients; dividing by X - a gives X - conj(a) back
    using P = polynomial<CQ>;
    using X_minus_a = P::sub_t<P::X, P::val<qa>>;
    using X_minus_conj_a = P::sub_t<P::X, P::val<CQ::conj_t<qa>>>;
    static_assert(std::is_same_v<P::div_t<P::mul_t<X_minus_a, X_minus_conj_a>, X_minus_a>, X_minus_conj_a>);
    // a X + 1 at i: -2 + 2^64 i
    EXPECT_EQ((P::val<qa, CQ::one>::eval(std::complex<double>(0, 1))), std::complex<double>(-2, 0x1p64));
}

TEST(complex, get_and_eval) {
    using cd = std::complex<double>;
    static_assert(q<1, 2, -3, 4>::get<cd>() == cd(0.5, -0.75));
    static_assert(q<1, 2, -3, 4>::get<std::complex<float>>() == std::complex<float>(0.5F, -0.75F));
    static_assert(z<3, -2>::get<cd>() == cd(3, -2) && f<5, 2, 4>::get<cd>() == cd(2, 4));
    // a complex value has no real conversion
    static_assert(!has_get<Qi::i, double> && has_get<Qi::i, cd>);
    // at run time: the multiplication of std::complex is not a constant expression with clang and libstdc++
    // (X - i)(X + i) = X² + 1 at 2 + i: (2 + i)² + 1 = 4 + 4i
    using P = polynomial<Qi>;
    using X2_plus_1 = P::mul_t<P::sub_t<P::X, P::val<Qi::i>>, P::add_t<P::X, P::val<Qi::i>>>;
    EXPECT_EQ(X2_plus_1::eval(cd(2, 1)), cd(4, 4));
    // (1/2 + i/4) X + (-1 + 2i) at -1 + 3i: -1/2 - 3/4 + (3/2 - 1/4) i - 1 + 2i = -9/4 + 13/4 i
    using L = P::val<q<1, 2, 1, 4>, q<-1, 1, 2, 1>>;
    EXPECT_EQ(L::eval(cd(-1, 3)), cd(-2.25, 3.25));
    const cd x(0.6, 0.8);
    EXPECT_EQ(L::eval(x), (q<1, 2, 1, 4>::get<cd>() * x + cd(-1, 2)));
}

TEST(complex, to_string) {
    EXPECT_EQ((z<3, -2>::to_string()), "3 - 2 i");
    EXPECT_EQ((z<3, 2>::to_string()), "3 + 2 i");
    EXPECT_EQ((z<3, 1>::to_string()), "3 + i");
    EXPECT_EQ((z<3, -1>::to_string()), "3 - i");
    EXPECT_EQ((z<0, 1>::to_string()), "i");
    EXPECT_EQ((z<0, -1>::to_string()), "-i");
    EXPECT_EQ((z<0, 5>::to_string()), "5 i");
    EXPECT_EQ((z<-4, 0>::to_string()), "-4");
    EXPECT_EQ((z<0, 0>::to_string()), "0");
    EXPECT_EQ((f<5, 2, 4>::to_string()), "2 + 4 i");
}
