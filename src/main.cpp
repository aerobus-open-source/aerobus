#include <cstdio>
#include <typeinfo>
#include <array>
#include <cmath>
#include <cstdlib>

#include "lib.h"


// conformance : https://godbolt.org/z/5chqjYEEs

using namespace aerobus;

int test_type_at() {
	if (!std::is_same<internal::type_at_t<0, float, int, long>, float>::value) {
		return 1;
	}
	if (!std::is_same<internal::type_at_t<1, float, int, long>, int>::value) {
		return 1;
	}
	if (!std::is_same<internal::type_at_t<2, float, int, long>, long>::value) {
		return 1;
	}

	return 0;
}

int test_poly_simplify() {
	using poly1 = polynomial<i32>::val<i32::val<0>, i32::val<1>, i32::val<2>>;
	using simplified1 = polynomial<i32>::simplify_t<poly1>;
	using expected1 = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
	if (!std::is_same<expected1, simplified1>::value) {
		return 1;
	}

	using poly2 = polynomial<i32>::val<i32::val<12>>;
	using simplified2 = polynomial<i32>::simplify_t<poly2>;
	if (!std::is_same<poly2, simplified2>::value) {
		return 1;
	}

	return 0;
}

int test_poly_eval() {
	// 1 + 2x + 3x^2
	using poly = polynomial<i32>::val<i32::val<3>, i32::val<2>, i32::val<1>>;
	constexpr int v = poly::eval(1);
	if (v != 6) {
		return 1;
	}
	constexpr float vv = poly::eval(1.0f);
	if (vv != 6.0f) {
		return 1;
	}

	// 1/2 + 3x/2
	using polyf = polynomial<Q32>::val<Q32::val<i32::val<3>, i32::val<2>>, Q32::val<i32::val<1>, i32::val<2>>>;
	constexpr float vvv = polyf::eval(1.0f);
	if (vvv != 2.0f) {
		return 1;
	}
	constexpr double vvvv = polyf::eval(-1.0);
	if (vvvv != -1.0) {
		return 1;
	}

	return 0;
}

int test_fraction_field_eval() {
	using half = Q32::val<i32::one, i32::val<2>>;
	constexpr float x = half::eval(2.0f);
	if (x != 0.5f) {
		return 1;
	}
	using thirdhalf = Q32::val<i32::val<3>, i32::val<2>>;
	constexpr float y = thirdhalf::eval(1.0f);
	if (y != 1.5f) {
		return 1;
	}

	// 3/2 + x / 2
	using polyA = polynomial<Q32>::val<half, thirdhalf>;
	constexpr float a = polyA::eval(2.0f);
	if (a != 2.5F) {
		return 1;
	}
	// 1/2 + x
	using polyB = polynomial<Q32>::val<Q32::one, half>;
	using F = FPQ32::val<polyA, polyB>;
	constexpr float z = F::eval(2.0f);
	if (z != 1.0f) {
		return 1;
	}
	constexpr float zz = F::eval(-1.0f);
	if (zz != -2.0f) {
		return 1;
	}


	return 0;
}

int test_coeff_at() {
	// 1 + 2x + 3x^2
	using poly = polynomial<i32>::val<i32::val<3>, i32::val<2>, i32::val<1>>;
	using at0 = poly::coeff_at_t<0>;
	using at1 = poly::coeff_at_t<1>;
	using at2 = poly::coeff_at_t<2>;
	using at3 = poly::coeff_at_t<3>;

	if(!std::is_same<at0, i32::val<1>>::value) {
		return 1;
	}

	if (!std::is_same<at1, i32::val<2>>::value) {
		return 1;
	}

	if (!std::is_same<at2, i32::val<3>>::value) {
		return 1;
	}

	if (!std::is_same<at3, i32::val<0>>::value) {
		return 1;
	}

	return 0;
}

template<typename... coeffs>
using IX = polynomial<i32>::val<coeffs...>;
template<int32_t x>
using Int = i32::val<x>;

template<typename P1, typename P2>
using add_ix = polynomial<i32>::add_t< P1, P2>;
template<typename P1, typename P2>
using sub_ix = polynomial<i32>::sub_t<P1, P2>;

int test_poly_add() {
	{
		// 1 + x
		using P1 = IX<Int<1>, Int<1>>;
		// 1 + x + x^2
		using P2 = IX<Int<1>, Int<1>, Int<1>>;
		// 2 + 2x + x^2
		using A = polynomial<i32>::add_t<P1, P2>;

		if (A::coeff_at_t<0>::v != 2) {
			return 1;
		}
		if (A::coeff_at_t<1>::v != 2) {
			return 1;
		}
		if (A::coeff_at_t<2>::v != 1) {
			return 1;
		}
		if (A::degree != 2) {
			return 1;
		}
	}

	{
		// 1 + x - x^2
		using P1 = polynomial<i32>::val<i32::val<-1>, i32::val<1>, i32::val<1>>;
		// 1 + x + x^2
		using P2 = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<1>>;
		// 2 + 2x
		using A = polynomial<i32>::add_t<P1, P2>;
		if (A::coeff_at_t<0>::v != 2) {
			return 1;
		}
		if (A::coeff_at_t<1>::v != 2) {
			return 1;
		}
		if (A::coeff_at_t<2>::v != 0) {
			return 1;
		}
		if (A::degree != 1) {
			return 1;
		}
	}

	return 0;
}

int test_poly_sub() {
	{
		// 1 + x
		using P1 = IX<Int<1>, Int<1>>;
		// 1 + x + x^2
		using P2 = IX<Int<1>, Int<1>, Int<1>>;
		// -x2
		using A = sub_ix<P1, P2>;

		if (A::coeff_at_t<0>::v != 0) {
			return 1;
		}
		if (A::coeff_at_t<1>::v != 0) {
			return 1;
		}
		if (A::coeff_at_t<2>::v != -1) {
			return 1;
		}
		if (A::degree != 2) {
			return 1;
		}
	}

	{
		// 1 + x + x^2
		using P1 = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<1>>;
		// 1 + x + x^2
		using P2 = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<1>>;
		// 0
		using A = polynomial<i32>::sub_t<P2, P1>;
		if (A::coeff_at_t<0>::v != 0) {
			return 1;
		}
		if (A::coeff_at_t<1>::v != 0) {
			return 1;
		}
		if (A::coeff_at_t<2>::v != 0) {
			return 1;
		}
		if (A::degree != 0) {
			return 1;
		}
	}

	return 0;
}

int test_poly_eq() {
	{
		using A = polynomial<i32>::val<i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>>;
		if (!polynomial<i32>::eq_t<A, B>::value) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<2>>;
		if (polynomial<i32>::eq_t<A, B>::value) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		if (!polynomial<i32>::eq_t<A, B>::value) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
		if (polynomial<i32>::eq_t<A, B>::value) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<2>>;
		if (polynomial<i32>::eq_t<A, B>::value) {
			return 1;
		}
	}

	return 0;
}

int test_gcd() {
	using A = i32::gcd_t<i32::val<12>, i32::val<6>>;
	if (A::v != 6) {
		return 1;
	}
	using B = i32::gcd_t<i32::val<12>, i32::val<6>>;
	if (B::v != 6) {
		return 1;
	}

	using C = i32::gcd_t<i32::val<5>, i32::val<3>>;
	if (C::v != 1) {
		return 1;
	}

	return 0;
}

int test_poly_mul() {
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<-1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using mul = polynomial<i32>::mul_t<A, B>;
		using expected = polynomial<i32>::val<i32::val<1>, i32::zero, i32::val<-1>>;
		if (!std::is_same<expected, mul>::value) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using mul = polynomial<i32>::mul_t<A, B>;
		using expected = polynomial<i32>::val<i32::val<1>, i32::val<2>, i32::val<1>>;
		if (!std::is_same<expected, mul>::value) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<2>>;
		using mul = polynomial<i32>::mul_t<A, B>;
		using expected = polynomial<i32>::val<i32::val<2>, i32::val<2>>;
		if (!std::is_same<expected, mul>::value) {
			return 1;
		}
	}

	return 0;
}

int test_monomial() {
	{
		//2x^3 + 0 + 0 + 0
		using A = polynomial<i32>::monomial_t<i32::val<2>, 3>;
		if (A::coeff_at_t<3>::v != 2) {
			return 1;
		}
		if (A::coeff_at_t<2>::v != 0) {
			return 1;
		}
		if (A::coeff_at_t<1>::v != 0) {
			return 1;
		}
		if (A::coeff_at_t<0>::v != 0) {
			return 1;
		}
	}

	return 0;
}

int test_poly_div() {
	// divisibility
	{
		// x^2 - 1
		using A = polynomial<i32>::val<i32::val<1>, i32::val<0>, i32::val<-1>>;
		// x - 1
		using B = polynomial<i32>::val<i32::val<1>, i32::val<-1>>;
		// x + 1
		using Q = polynomial<i32>::div_t<A, B>;
		using R = polynomial<i32>::mod_t<A, B>;
		if (!R::is_zero_t::value) {
			return 1;
		}
		if (Q::degree != 1) {
			return 1;
		}
		if (Q::coeff_at_t<0>::v != 1) {
			return 1;
		}
		if (Q::coeff_at_t<1>::v != 1) {
			return 1;
		}
	}
	// divide by constant
	{
		using A = polynomial<i32>::val<i32::val<2>, i32::val<2>, i32::val<2>>;
		using B = polynomial<i32>::val<i32::val<2>>;
		using C = polynomial<i32>::div_t<A, B>;
		using R = polynomial<i32>::mod_t<A, B>;
		if (!R::is_zero_t::value) {
			return 1;
		}
		if (C::degree != 2) {
			return 1;
		}
		if (C::coeff_at_t<0>::v != 1) {
			return 1;
		}
		if (C::coeff_at_t<1>::v != 1) {
			return 1;
		}
		if (C::coeff_at_t<2>::v != 1) {
			return 1;
		}
	}
	// no divisibility
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using C = polynomial<i32>::div_t<A, B>;
		using R = polynomial<i32>::mod_t<A, B>;
		if (C::degree != 1) {
			return 1;
		}
		if (C::coeff_at_t<0>::v != 0) {
			return 1;
		}
		if (C::coeff_at_t<1>::v != 1) {
			return 1;
		}
		if (R::degree != 0) {
			return 1;
		}
		if (R::aN::v != 1) {
			return 1;
		}
	}
	// float divisibility
	{
		// x2 -1
		using A = polynomial<Q32>::val<Q32::one, Q32::zero, Q32::val<i32::val<-1>, i32::val<1>>>;
		// 2x + 2
		using B = polynomial<Q32>::val<Q32::val<i32::val<2>, i32::one>, Q32::val<i32::val<2>, i32::one>>;
		using Q = polynomial<Q32>::div_t<A, B>;

		if (Q::degree != 1) {
			return 1;
		}
		if (Q::coeff_at_t<0>::x::v != -1) {
			return 1;
		}
		if (Q::coeff_at_t<0>::y::v != 2) {
			return 1;
		}
		if (Q::coeff_at_t<1>::x::v != 1) {
			return 1;
		}
		if (Q::coeff_at_t<1>::y::v != 2) {
			return 1;
		}
	}
	// float divisibility
	{
		// x2 -1
		using A = polynomial<Q32>::val<Q32::one, Q32::one>;
		// 2x + 2
		using B = polynomial<Q32>::val<Q32::val<i32::val<2>, i32::one>, Q32::val<i32::val<2>, i32::one>>;
		using Q = polynomial<Q32>::div_t<A, B>;
		if (Q::degree != 0) {
			return 1;
		}
		if (Q::coeff_at_t<0>::x::v != 1) {
			return 1;
		}
		if (Q::coeff_at_t<0>::y::v != 2) {
			return 1;
		}
	}
	return 0;
}

int test_poly_gcd() {
	{
		// (x+1)*(x+1)
		using A = polynomial<Q32>::val<Q32::one, Q32::val<i32::val<2>, i32::val<1>>, Q32::one>;
		// (x+1)
		using B = polynomial<Q32>::val<Q32::one, Q32::one>;
		using G = internal::gcd<polynomial<Q32>>::type<A, B>;
		if (G::degree != 1) {
			return 1;
		}
		if (G::coeff_at_t<0>::x::v != 1 || G::coeff_at_t<0>::y::v != 1) {
			return 1;
		}
		if (G::coeff_at_t<1>::x::v != 1 || G::coeff_at_t<1>::y::v != 1) {
			return 1;
		}
	}
	{
		// (x+1)*(x+1)
		using A = polynomial<Q32>::val<Q32::one, Q32::val<i32::val<2>, i32::val<1>>, Q32::one>;
		// (x+1)*(x-1) :: x^2 - 1
		using B = polynomial<Q32>::val<Q32::one, Q32::zero, Q32::val<i32::val<-1>, i32::val<1>>>;
		// x + 1
		using G = polynomial<Q32>::gcd_t<A, B>;
		if (G::degree != 1) {
			return 1;
		}
		if (G::coeff_at_t<0>::x::v != 1 || G::coeff_at_t<0>::y::v != 1) {
			return 1;
		}
		if (G::coeff_at_t<1>::x::v != 1 || G::coeff_at_t<1>::y::v != 1) {
			return 1;
		}
	}

	return 0;
}

int test_poly_to_string() {
	using P32 = polynomial<Q32>;
	using A = FPQ32::val<P32::val<Q32::one, Q32::one>, P32::val<Q32::val<i32::val<2>, i32::val<1>>, Q32::one>>;
	auto rep = A::to_string();
	const char* expected = "(x + 1) / (2 x + 1)";
	if (strcmp(rep.c_str(),expected) != 0) {
		printf("expected %s got %s\n", "(x + 1) / (2 x + 1)", rep.c_str());
		return 1;
	}

	return 0;
}

int test_add_q32() {
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<4>>;

		using c = Q32::add_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 3 || y != 4) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<3>>;

		using c = Q32::add_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 5 || y != 6) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<-1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<2>>;

		using c = Q32::add_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 0 || y != 1) {
			return 1;
		}
	}

	return 0;
}

int test_sub_q32() {
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<4>>;

		using c = Q32::sub_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 4) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<3>>;

		using c = Q32::sub_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 6) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<2>>;

		using c = Q32::sub_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 0 || y != 1) {
			return 1;
		}
	}

	return 0;
}

int test_mul_q32() {
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<4>>;

		using c = Q32::mul_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 8) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<3>>;

		using c = Q32::mul_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 6) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<0>, i32::val<2>>;

		using c = Q32::mul_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 0 || y != 1) {
			return 1;
		}
	}

	return 0;
}

int test_div_q32() {
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<4>>;

		using c = Q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 2 || y != 1) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<1>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<3>>;

		using c = Q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 3 || y != 2) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<0>, i32::val<2>>;
		using b = Q32::val<i32::val<1>, i32::val<2>>;

		using c = Q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 0 || y != 1) {
			return 1;
		}
	}
	{
		using a = Q32::val<i32::val<1>, i32::val<1>>;
		using b = Q32::val<i32::val<2>, i32::val<1>>;

		using c = Q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 2) {
			return 1;
		}
	}

	return 0;
}

int test_simplify_q32() {
	using A = Q32::val<i32::val<2>, i32::val<2>>;
	using B = Q32::val<i32::val<1>, i32::val<1>>;
	using C = Q32::val<i32::val<-1>, i32::val<-1>>;
	using D = Q32::val<i32::val<1>, i32::val<-2>>;
	if (!Q32::eq_t<Q32::simplify_t<A>, Q32::one>::value) {
		return 1;
	}
	if (!Q32::eq_t<Q32::simplify_t<B>, Q32::one>::value) {
		return 1;
	}
	if (!Q32::eq_t<Q32::simplify_t<C>, Q32::one>::value) {
		return 1;
	}
	if (!Q32::eq_t<Q32::simplify_t<D>, Q32::val<i32::val<-1>, i32::val<2>>>::value) {
		return 1;
	}
	if(!Q32::eq_t<Q32::simplify_t<Q32::sub_t<Q32::zero, Q32::zero>>, Q32::zero>::value) {
		return 1;
	}

	return 0;
}

int test_eq_q32() {
	using A = Q32::val<i32::val<2>, i32::val<2>>;
	using B = Q32::val<i32::val<1>, i32::val<1>>;
	using C = Q32::val<i32::val<-1>, i32::val<-1>>;
	if (!Q32::eq_t<A, B>::value) {
		return 1;
	}
	if (!Q32::eq_t<A, C>::value) {
		return 1;
	}
	return 0;
}

int test_fraction_field_of_fraction_field () {
	using qq32 = FractionField<Q32>;
	if (!std::is_same<Q32, qq32>::value) {
		return 1;
	}

	return 0;
}

int test_factorial() {
	constexpr float x = factorial<i32, 3>::value;
	if (x != 6.0f) {
		return 1;
	}

	constexpr size_t y = factorial<i32, 0>::value;
	if (y != 1) {
		return 1;
	}

	return 0;
}

int test_combination() {
	constexpr int x = combination<i32, 2, 4>::value;
	if (x != 6) {
		return 1;
	}
	constexpr int y = combination<i32, 0, 4>::value;
	if (y != 1) {
		return 1;
	}
	constexpr int z = combination<i32, 1, 4>::value;
	if (z != 4) {
		return 1;
	}

	constexpr int zz = combination<i32, 3, 4>::value;
	if (zz != 4) {
		return 1;
	}
	return 0;
}

int test_bernouilli() {
	constexpr float b0 = bernouilli<i32, 0>::template value<float>;
	if (b0 != 1.0f) {
		return 1;
	}
	constexpr float b1 = bernouilli<i32, 1>::template value<float>;
	if (b1 != -0.5f) {
		return 1;
	}
	using B2 = bernouilli<i32, 2>::type;
	if (B2::x::v != 1 || B2::y::v != 6) {
		return 1;
	}
	constexpr double b3 = bernouilli<i32, 3>::template value<double>;
	if (b3 != 0.0) {
		return 1;
	}

	return 0;
}

int test_zpz () {
	using Z2Z = zpz<2>;
	if (Z2Z::template add_t<typename Z2Z::val<1>, typename Z2Z::val<1>>::v != 0) {
		return 1;
	}
	if (!Z2Z::is_field) {
		return 1;
	}

	using Z4Z = zpz<4>;
	if (Z4Z::template add_t<typename Z4Z::val<4>, typename Z4Z::val<12>>::v != 0) {
		return 1;
	}
	if (Z4Z::template add_t<typename Z4Z::val<5>, typename Z4Z::val<13>>::v == 0) {
		return 1;
	}
	if (Z4Z::is_field) {
		return 1;
	}

	// gcd
	if (Z2Z::template gcd_t<typename Z2Z::val<1>, typename Z2Z::val<1>>::v != 1) {
		return 1;
	}
	using Z5Z = zpz<5>;
	if (Z5Z::template gcd_t<typename Z5Z::val<2>, typename Z5Z::val<3>>::v != 1) {
		return 1;
	}
	if (Z5Z::template gcd_t<typename Z5Z::val<2>, typename Z5Z::val<4>>::v != 2) {
		return 1;
	}

	return 0;
}

int test_is_prime() {
	if (is_prime<1>::value) {
		return 1;
	}
	if (!is_prime<2>::value) {
		return 1;
	}
	if (!is_prime<3>::value) {
		return 1;
	}
	if (is_prime<4>::value) {
		return 1;
	}
	if (!is_prime<5>::value) {
		return 1;
	}
	if (is_prime<6>::value) {
		return 1;
	}
	if (!is_prime<7>::value) {
		return 1;
	}
	if (is_prime<8>::value) {
		return 1;
	}
	if (is_prime<9>::value) {
		return 1;
	}
	if (is_prime<10>::value) {
		return 1;
	}
	if (!is_prime<31>::value) {
		return 1;
	}
	if (is_prime<100>::value) {
		return 1;
	}

	return 0;
}

int test_exp() {
	using E = aerobus::exp<i32, 12>;
	constexpr float e0 = E::eval(0.0F);
	constexpr float e01 = E::eval(0.1F);
	if (e0 != 1.0F) {
		return 1;
	}

	if (std::abs(std::exp(0.1f) - e01) > 1E-7F) {
		return 1;
	}

	return 0;
}

int test_alternate() {
	constexpr int a0 = alternate<i32, 0>::value;
	if (a0 != 1) {
		return 1;
	}
	constexpr int a1 = alternate<i32, 1>::value;
	if (a1 != -1) {
		return 1;
	}
	constexpr int a2 = alternate<i32, 2>::value;
	if (a2 != 1) {
		return 1;
	}

	return 0;
}

int test_type_list () {
	using A = type_list<int, float>;
	if (A::length != 2) {
		return 1;
	}
	if (typeid(A::at<0>) != typeid(int)) {
		return 1;
	}
	if (typeid(A::at<1>) != typeid(float)) {
		return 1;
	}

	using B = A::push_back<double>;
	if (B::length != 3) {
		return 1;
	}
	if(typeid(B::at<2>) != typeid(double)) {
		return 1;
	}

	using C = B::pop_front;
	if(typeid(C::type) != typeid(int)) {
		return 1;
	}
	if(C::tail::length != 2) {
		return 1;
	}
	if (typeid(C::tail::at<0>) != typeid(float)) {
		return 1;
	}
	using D = C::tail::split<1>;
	if(D::head::length != 1 || D::tail::length != 1) {
		return 1;
	}
	if(typeid(D::head::at<0>) != typeid(float) || typeid(D::tail::at<0>) != typeid(double)) {
		return 1;
	}

	return 0;
}

int test_concept_ring() {
	static_assert(aerobus::IsRing<i32>);
	static_assert(aerobus::IsRing<i64>);
	static_assert(aerobus::IsRing<zpz<3>>);
	static_assert(aerobus::IsRing<Q32>);
	static_assert(aerobus::IsField<Q32>);
	static_assert(aerobus::IsRing<polynomial<i32>>);
	static_assert(aerobus::IsField<FractionField<polynomial<i32>>>);
	return 0;
}

#define RUN_TEST(test_name) if (test_name() != 0) { printf("%s failed\n", #test_name); return 1; }

int main(int argc, char* argv[]) {
	RUN_TEST(test_concept_ring)
	RUN_TEST(test_type_list)
	RUN_TEST(test_type_at) 
	RUN_TEST(test_poly_simplify)
	RUN_TEST(test_coeff_at)
	RUN_TEST(test_poly_add)
	RUN_TEST(test_poly_sub)
	RUN_TEST(test_poly_eq)
	RUN_TEST(test_gcd)
	RUN_TEST(test_poly_mul)
	RUN_TEST(test_poly_to_string)
	RUN_TEST(test_monomial)
	RUN_TEST(test_poly_gcd)
	RUN_TEST(test_poly_eval)
	RUN_TEST(test_add_q32)
	RUN_TEST(test_fraction_field_eval)
	RUN_TEST(test_sub_q32)
	RUN_TEST(test_mul_q32)
	RUN_TEST(test_div_q32)
	RUN_TEST(test_eq_q32)
	RUN_TEST(test_fraction_field_of_fraction_field)
	RUN_TEST(test_factorial)
	RUN_TEST(test_combination)
	RUN_TEST(test_bernouilli)
	RUN_TEST(test_alternate)
	RUN_TEST(test_exp)
	RUN_TEST(test_is_prime)
	RUN_TEST(test_zpz)
	printf("ALL TESTS OK\n");
	return 0;
}