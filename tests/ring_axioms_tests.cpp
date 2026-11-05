#include <gtest/gtest.h>
#include <type_traits>
#include <utility>

#define AEROBUS_CONWAY_IMPORTS
#include "aerobus.h"

using namespace aerobus;  // NOLINT

namespace {
// a ring and a few of its values
template<typename Ring, typename... Values>
struct samples {
    using ring = Ring;
    using values = type_list<Values...>;
    static constexpr size_t count = sizeof...(Values);
};

template<typename S, size_t i>
using sample = typename S::values::template at<i>;

// the axioms for one triple of values; equality is value-based (eq_v)
template<typename S, size_t i, size_t j, size_t k>
constexpr bool axioms_hold() {
    using R = typename S::ring;
    using a = sample<S, i>;
    using b = sample<S, j>;
    using c = sample<S, k>;
    using add = typename R::template add_t<a, b>;
    using mul = typename R::template mul_t<a, b>;
    bool ok = R::template eq_v<add, typename R::template add_t<b, a>>
        && R::template eq_v<mul, typename R::template mul_t<b, a>>
        && R::template eq_v<typename R::template add_t<add, c>,
                            typename R::template add_t<a, typename R::template add_t<b, c>>>
        && R::template eq_v<typename R::template mul_t<mul, c>,
                            typename R::template mul_t<a, typename R::template mul_t<b, c>>>
        && R::template eq_v<typename R::template mul_t<a, typename R::template add_t<b, c>>,
                            typename R::template add_t<mul, typename R::template mul_t<a, c>>>
        && R::template eq_v<typename R::template add_t<a, typename R::zero>, a>
        && R::template eq_v<typename R::template mul_t<a, typename R::one>, a>
        && R::template eq_v<typename R::template add_t<typename R::template sub_t<typename R::zero, a>, a>,
                            typename R::zero>;
    if constexpr (R::is_field) {
        // in a field, division by a nonzero value is exact
        if constexpr (!b::is_zero_t::value) {
            ok = ok && R::template eq_v<typename R::template mul_t<typename R::template div_t<a, b>, b>, a>;
        }
    }
    return ok;
}

template<typename S, size_t... I>
constexpr bool all_triples(std::index_sequence<I...>) {
    constexpr size_t n = S::count;
    return (axioms_hold<S, I / (n * n), (I / n) % n, I % n>() && ...);
}

template<typename S>
constexpr bool ring_axioms() {
    return all_triples<S>(std::make_index_sequence<S::count * S::count * S::count>{});
}

using F9 = GaloisField<3, 2>;
using phi = F9::inject_ring_t<polynomial<zpz<3>>::X>;
using Z17 = Quotient<i64, i64::val<17>>;

template<typename S>
class ring_axioms_test : public ::testing::Test {};

using all_samples = ::testing::Types<
    samples<i32, i32::val<-7>, i32::zero, i32::val<12>>,
    samples<i64, i64::val<-3>, i64::one, i64::val<100000>>,
    samples<zpz<7>, zpz<7>::val<-1>, zpz<7>::val<3>, zpz<7>::zero>,
    samples<zpz<8>, zpz<8>::val<2>, zpz<8>::val<4>, zpz<8>::val<7>>,
    samples<q32, make_q32_t<-1, 2>, make_q32_t<2, 3>, q32::zero>,
    samples<pi64, make_int_polynomial_t<i64, 1, -1>, make_int_polynomial_t<i64, 2, 0, 3>, pi64::one>,
    samples<pq64, make_frac_polynomial_t<i64, 1, 2>, pq64::val<make_q64_t<1, 3>>, pq64::X>,
    samples<Z17, Z17::inject_constant_t<-1>, Z17::inject_constant_t<5>, Z17::inject_constant_t<9>>,
    samples<F9, phi, F9::mul_t<phi, phi>, F9::inject_constant_t<2>>,
    samples<dyadic, dyadic::from_hex_float_t<"-0x3p-2">, dyadic::from_hex_float_t<"0x5p+40">,
            dyadic::from_hex_float_t<"0x1p-90">>,
    samples<polynomial<dyadic>, polynomial<dyadic>::val<dyadic::one, dyadic::from_hex_float_t<"0x1p-3">>,
            polynomial<dyadic>::X, polynomial<dyadic>::one>>;
}  // namespace

TYPED_TEST_SUITE(ring_axioms_test, all_samples);

TYPED_TEST(ring_axioms_test, hold_on_samples) {
    EXPECT_TRUE(ring_axioms<TypeParam>());
}
