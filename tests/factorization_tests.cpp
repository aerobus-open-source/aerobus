#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#define AEROBUS_CONWAY_IMPORTS
#include "aerobus.h"

// Expected factorizations over F_p were computed with SymPy (factor_list(..., modulus=p)). Factors are
// monic and sorted by the order of polynomials: degree first, then coefficients from the highest.

namespace {
using aerobus::factor;
using aerobus::factorization;
using aerobus::type_list;

template<int32_t p, int32_t... coeffs>
using P = aerobus::make_int_polynomial_t<aerobus::zpz<p>, coeffs...>;

template<int32_t p>
using X = typename aerobus::polynomial<aerobus::zpz<p>>::X;

template<int64_t... coeffs>
using PQ = aerobus::embed_int_poly_in_fractions_t<aerobus::make_int_polynomial_t<aerobus::i64, coeffs...>>;

// x^n - 1 with coefficients in F
template<typename F, size_t n>
using xn_minus_one = typename aerobus::polynomial<F>::template sub_t<
    typename aerobus::polynomial<F>::template monomial_t<typename F::one, n>,
    typename aerobus::polynomial<F>::one>;

template<typename Fact>
constexpr bool all_factors_irreducible() {
    return []<typename... Fs>(type_list<Fs...>) {
        return (aerobus::is_irreducible_v<typename Fs::type> && ...);
    }(typename Fact::factors{});
}

template<typename Fact, size_t degree>
constexpr size_t count_of_degree() {
    return []<typename... Fs>(type_list<Fs...>) {
        return ((Fs::type::degree == degree ? size_t{1} : size_t{0}) + ... + size_t{0});
    }(typename Fact::factors{});
}
}  // namespace

TEST(finite_fields, concept_and_traits) {
    static_assert(aerobus::IsFiniteField<aerobus::zpz<7>>);
    static_assert(!aerobus::IsFiniteField<aerobus::zpz<8>>);
    static_assert(!aerobus::IsFiniteField<aerobus::i64>);
    static_assert(!aerobus::IsFiniteField<aerobus::q64>);
    static_assert(aerobus::IsFiniteField<aerobus::GaloisField<2, 8>>);
    static_assert(aerobus::IsFiniteField<aerobus::Quotient<aerobus::i64, aerobus::i64::val<17>>>);
    static_assert(aerobus::characteristic_v<aerobus::GaloisField<2, 8>> == 2);
    static_assert(aerobus::field_order_v<aerobus::GaloisField<2, 8>> == 256);
    static_assert(aerobus::field_order_v<aerobus::GaloisField<3, 2>> == 9);
    static_assert(aerobus::field_order_v<aerobus::zpz<13>> == 13);
    static_assert(aerobus::field_order_v<aerobus::Quotient<aerobus::i64, aerobus::i64::val<17>>> == 17);
}

TEST(finite_fields, elements_are_enumerated) {
    using F9 = aerobus::GaloisField<3, 2>;
    // element i has the base-p digits of i as coefficients, lowest first: 5 = 2 + 1 * 3 -> x + 2
    static_assert(std::is_same_v<aerobus::field_element_t<F9, 5>, F9::inject_ring_t<P<3, 1, 2>>>);
    static_assert(std::is_same_v<aerobus::field_element_t<F9, 0>, F9::zero>);
    static_assert(std::is_same_v<aerobus::field_element_t<aerobus::zpz<7>, 4>, aerobus::zpz<7>::val<4>>);
}

TEST(polynomials, monic) {
    static_assert(std::is_same_v<aerobus::monic_t<P<5, 3, 0, 3>>, P<5, 1, 0, 1>>);
    static_assert(std::is_same_v<aerobus::monic_t<PQ<2, 1>>, aerobus::pq64::val<aerobus::q64::one,
                                                                                aerobus::make_q64_t<1, 2>>>);
}

TEST(irreducibility, prime_fields) {
    static_assert(aerobus::is_irreducible_v<P<2, 1, 1>>);               // x + 1
    static_assert(aerobus::is_irreducible_v<P<2, 1, 1, 1>>);            // x^2 + x + 1
    static_assert(!aerobus::is_irreducible_v<P<2, 1, 0, 1>>);           // x^2 + 1 = (x + 1)^2
    static_assert(aerobus::is_irreducible_v<P<3, 1, 0, 1>>);            // x^2 + 1 over F3
    static_assert(!aerobus::is_irreducible_v<P<3, 1, 0, 0, 0, 1>>);     // x^4 + 1 = (x^2 + x + 2)(x^2 + 2x + 2)
    static_assert(!aerobus::is_irreducible_v<P<5, 1, 0, 0, 0, 0, 4>>);  // x^5 - 1 = (x + 4)^5
    static_assert(!aerobus::is_irreducible_v<P<2, 1>>);                 // constants are units
}

TEST(irreducibility, conway_polynomials) {
    static_assert(aerobus::is_irreducible_v<aerobus::ConwayPolynomial<2, 8>::type>);
    static_assert(aerobus::is_irreducible_v<aerobus::ConwayPolynomial<3, 5>::type>);
    static_assert(aerobus::is_irreducible_v<aerobus::ConwayPolynomial<5, 4>::type>);
    static_assert(aerobus::is_irreducible_v<aerobus::ConwayPolynomial<7, 3>::type>);
    static_assert(aerobus::is_irreducible_v<aerobus::ConwayPolynomial<2, 12>::type>);
}

TEST(irreducibility, extension_fields) {
    using F9 = aerobus::GaloisField<3, 2>;
    using PF9 = aerobus::polynomial<F9>;
    using phi = F9::inject_ring_t<X<3>>;  // a generator of F9*, hence not a square
    static_assert(aerobus::is_irreducible_v<PF9::val<F9::one, F9::zero, F9::sub_t<F9::zero, phi>>>);      // x^2 - phi
    static_assert(!aerobus::is_irreducible_v<PF9::val<F9::one, F9::zero, F9::sub_t<F9::zero, F9::one>>>);  // x^2 - 1
}

TEST(squarefree, prime_fields) {
    // (x + 1)^2 (x^2 + x + 1)
    static_assert(std::is_same_v<aerobus::squarefree_t<P<2, 1, 1, 0, 1, 1>>,
                                 factorization<aerobus::zpz<2>::one,
                                               type_list<factor<P<2, 1, 1>, 2>, factor<P<2, 1, 1, 1>, 1>>>>);
    // (x + 2)(x^2 + 1)^3 over F3: the cube is a polynomial in x^3, found through a p-th root
    static_assert(std::is_same_v<aerobus::squarefree_t<P<3, 1, 2, 0, 0, 0, 0, 1, 2>>,
                                 factorization<aerobus::zpz<3>::one,
                                               type_list<factor<P<3, 1, 2>, 1>, factor<P<3, 1, 0, 1>, 3>>>>);
    // x^5 - 1 = (x + 4)^5 over F5
    static_assert(std::is_same_v<aerobus::squarefree_t<P<5, 1, 0, 0, 0, 0, 4>>,
                                 factorization<aerobus::zpz<5>::one, type_list<factor<P<5, 1, 4>, 5>>>>);
    // already squarefree, not irreducible: x^4 + 1 over F3
    static_assert(std::is_same_v<aerobus::squarefree_t<P<3, 1, 0, 0, 0, 1>>,
                                 factorization<aerobus::zpz<3>::one, type_list<factor<P<3, 1, 0, 0, 0, 1>, 1>>>>);
}

TEST(squarefree, rationals) {
    // (x + 1)^2 (x - 2)
    static_assert(std::is_same_v<aerobus::squarefree_t<PQ<1, 0, -3, -2>>,
                                 factorization<aerobus::q64::one,
                                               type_list<factor<PQ<1, -2>, 1>, factor<PQ<1, 1>, 2>>>>);
    // the unit is the leading coefficient
    static_assert(std::is_same_v<typename aerobus::squarefree_t<PQ<2, 4, 2>>::unit, aerobus::make_q64_t<2, 1>>);
}

TEST(factorization, prime_fields) {
    // x^7 - 1 over F2
    static_assert(std::is_same_v<aerobus::factorization_t<xn_minus_one<aerobus::zpz<2>, 7>>,
                                 factorization<aerobus::zpz<2>::one,
                                               type_list<factor<P<2, 1, 1>, 1>, factor<P<2, 1, 0, 1, 1>, 1>,
                                                         factor<P<2, 1, 1, 0, 1>, 1>>>>);
    // x^15 - 1 over F2: the minimal polynomials of the BCH codes of length 15
    static_assert(std::is_same_v<aerobus::factorization_t<xn_minus_one<aerobus::zpz<2>, 15>>,
                                 factorization<aerobus::zpz<2>::one,
                                               type_list<factor<P<2, 1, 1>, 1>, factor<P<2, 1, 1, 1>, 1>,
                                                         factor<P<2, 1, 0, 0, 1, 1>, 1>, factor<P<2, 1, 1, 0, 0, 1>, 1>,
                                                         factor<P<2, 1, 1, 1, 1, 1>, 1>>>>);
    // x^8 - 1 over F3
    static_assert(std::is_same_v<aerobus::factorization_t<xn_minus_one<aerobus::zpz<3>, 8>>,
                                 factorization<aerobus::zpz<3>::one,
                                               type_list<factor<P<3, 1, 1>, 1>, factor<P<3, 1, 2>, 1>,
                                                         factor<P<3, 1, 0, 1>, 1>, factor<P<3, 1, 1, 2>, 1>,
                                                         factor<P<3, 1, 2, 2>, 1>>>>);
    // (x + 2)(x^2 + 1)^3 over F3: multiplicities are kept
    static_assert(std::is_same_v<aerobus::factorization_t<P<3, 1, 2, 0, 0, 0, 0, 1, 2>>,
                                 factorization<aerobus::zpz<3>::one,
                                               type_list<factor<P<3, 1, 2>, 1>, factor<P<3, 1, 0, 1>, 3>>>>);
    // 3x^4 + 3 over F5: the unit is the leading coefficient
    static_assert(std::is_same_v<aerobus::factorization_t<P<5, 3, 0, 0, 0, 3>>,
                                 factorization<aerobus::zpz<5>::val<3>,
                                               type_list<factor<P<5, 1, 0, 2>, 1>, factor<P<5, 1, 0, 3>, 1>>>>);
}

TEST(factorization, product_is_the_polynomial) {
    using A = P<3, 1, 2, 0, 0, 0, 0, 1, 2>;
    static_assert(std::is_same_v<typename aerobus::factorization_t<A>::type, A>);
    using B = P<5, 3, 0, 0, 0, 3>;
    static_assert(std::is_same_v<typename aerobus::factorization_t<B>::type, B>);
}

TEST(factorization, extension_fields) {
    // x^8 - 1 splits into linear factors over F9, whose multiplicative group has order 8
    using F9 = aerobus::GaloisField<3, 2>;
    using A = aerobus::factorization_t<xn_minus_one<F9, 8>>;
    static_assert(A::factors::length == 8 && count_of_degree<A, 1>() == 8);
    static_assert(std::is_same_v<typename A::type, xn_minus_one<F9, 8>>);
    // x^15 - 1 over F4: cyclotomic cosets of 4 modulo 15 give 3 linear and 6 quadratic factors
    using F4 = aerobus::GaloisField<2, 2>;
    using B = aerobus::factorization_t<xn_minus_one<F4, 15>>;
    static_assert(count_of_degree<B, 1>() == 3 && count_of_degree<B, 2>() == 6 && B::factors::length == 9);
    static_assert(all_factors_irreducible<B>());
    static_assert(std::is_same_v<typename B::type, xn_minus_one<F4, 15>>);
}
