#include <cstdio>
#include <typeinfo>
#include <array>
#include <cmath>
#include <cstdlib>

#include "lib.h"
#include "../imports/conwaypolynomials.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

// conformance : https://godbolt.org/z/qnfP99KWv

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
	using polyf = polynomial<q32>::val<q32::val<i32::val<3>, i32::val<2>>, q32::val<i32::val<1>, i32::val<2>>>;
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
	using half = q32::val<i32::one, i32::val<2>>;
	constexpr float x = half::eval(2.0f);
	if (x != 0.5f) {
		return 1;
	}
	using thirdhalf = q32::val<i32::val<3>, i32::val<2>>;
	constexpr float y = thirdhalf::eval(1.0f);
	if (y != 1.5f) {
		return 1;
	}

	// 3/2 + x / 2
	using polyA = polynomial<q32>::val<half, thirdhalf>;
	constexpr float a = polyA::eval(2.0f);
	if (a != 2.5F) {
		return 1;
	}
	// 1/2 + x
	using polyB = polynomial<q32>::val<q32::one, half>;
	using F = fpq32::val<polyA, polyB>;
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

int test_poly_derive() {
	{
		// 1 + x
		using P1 = IX<Int<1>, Int<1>>;
		using PP = polynomial<i32>::template derive_t<P1>;
		if(PP::degree != 0) {
			return 1;
		}
		if(PP::coeff_at_t<0>::v != 1) {
			return 1;
		}
	} 
	{
		// 1 + x + 3x�
		using P1 = IX<Int<3>, Int<1>, Int<1>>;
		using PP = polynomial<i32>::template derive_t<P1>;
		if (PP::degree != 1) {
			return 1;
		}
		if (PP::coeff_at_t<0>::v != 1) {
			return 1;
		}
		if (PP::coeff_at_t<1>::v != 6) {
			return 1;
		}
	} 
	using Z2Z = zpz<2>;
	{
		// in Z/2Z
		// 1 + x + 3x� -> 1
		using P1 = polynomial<Z2Z>::template val<Z2Z::template val<3>, Z2Z::template val<1>, Z2Z::template val<1>>;
		using PP = polynomial<Z2Z>::template derive_t<P1>;
		if (PP::degree != 0) {
			return 1;
		}
		if (PP::coeff_at_t<0>::v != 1) {
			return 1;
		}
	}
	{
		// in Z/2Z
		// x^3 + x^2 + x + 1 -> 1 + x^2 (because 3 == 1)
		using P1 = polynomial<Z2Z>::template val<Z2Z::template val<1>, Z2Z::template val<1>, Z2Z::template val<1>, Z2Z::template val<1>>;
		using PP = polynomial<Z2Z>::template derive_t<P1>;
		if (PP::degree != 2) {
			return 1;
		}
		if (PP::coeff_at_t<0>::v != 1) {
			return 1;
		}
		if (PP::coeff_at_t<1>::v != 0) {
			return 1;
		}
		if (PP::coeff_at_t<2>::v != 1) {
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
		if (!polynomial<i32>::eq_v<A, B>) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<2>>;
		if (polynomial<i32>::eq_v<A, B>) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		if (!polynomial<i32>::eq_v<A, B>) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
		if (polynomial<i32>::eq_v<A, B>) {
			return 1;
		}
	}
	{
		using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
		using B = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<2>>;
		if (polynomial<i32>::eq_v<A, B>) {
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
		if (!R::is_zero_v) {
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
		if (!R::is_zero_v) {
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
		using A = polynomial<q32>::val<q32::one, q32::zero, q32::val<i32::val<-1>, i32::val<1>>>;
		// 2x + 2
		using B = polynomial<q32>::val<q32::val<i32::val<2>, i32::one>, q32::val<i32::val<2>, i32::one>>;
		using Q = polynomial<q32>::div_t<A, B>;

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
		using A = polynomial<q32>::val<q32::one, q32::one>;
		// 2x + 2
		using B = polynomial<q32>::val<q32::val<i32::val<2>, i32::one>, q32::val<i32::val<2>, i32::one>>;
		using Q = polynomial<q32>::div_t<A, B>;
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
		using A = polynomial<q32>::val<q32::one, q32::val<i32::val<2>, i32::val<1>>, q32::one>;
		// (x+1)
		using B = polynomial<q32>::val<q32::one, q32::one>;
		using G = internal::gcd<polynomial<q32>>::type<A, B>;
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
		using A = polynomial<q32>::val<q32::one, q32::val<i32::val<2>, i32::val<1>>, q32::one>;
		// (x+1)*(x-1) :: x^2 - 1
		using B = polynomial<q32>::val<q32::one, q32::zero, q32::val<i32::val<-1>, i32::val<1>>>;
		// x + 1
		using G = polynomial<q32>::gcd_t<A, B>;
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
	using P32 = polynomial<q32>;
	using A = fpq32::val<P32::val<q32::one, q32::one>, P32::val<q32::val<i32::val<2>, i32::val<1>>, q32::one>>;
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
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<4>>;

		using c = q32::add_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 3 || y != 4) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<3>>;

		using c = q32::add_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 5 || y != 6) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<-1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<2>>;

		using c = q32::add_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 0 || y != 1) {
			return 1;
		}
	}

	return 0;
}

int test_is_zero_q32() {
	using a = q32::zero;
	if (!a::is_zero_v) {
		return 1;
	}
	using b = q32::one;
	if (b::is_zero_v) {
		return 1;
	}

	return 0;
}

int test_gt_q32() {
	{
		using a = q32::inject_constant_t<1>;
		using b = q32::zero;
		if (!q32::gt_v<a, b>) {
			return 1;
		}
	} 
	{
		using a = q32::zero;
		using b = q32::inject_constant_t<2>;
		if (q32::gt_v<a, b>) {
			return 1;
		}
	}
	{
		using a = q32::zero;
		using b = q32::zero;
		if (q32::gt_v<a, b>) {
			return 1;
		}
	}

	return 0;
}

int test_sub_q32() {
	{
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<4>>;

		using c = q32::sub_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 4) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<3>>;

		using c = q32::sub_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 6) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<2>>;

		using c = q32::sub_t<a, b>;

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
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<4>>;

		using c = q32::mul_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 8) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<3>>;

		using c = q32::mul_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 6) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<0>, i32::val<2>>;

		using c = q32::mul_t<a, b>;

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
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<4>>;

		using c = q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 2 || y != 1) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<1>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<3>>;

		using c = q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 3 || y != 2) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<0>, i32::val<2>>;
		using b = q32::val<i32::val<1>, i32::val<2>>;

		using c = q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 0 || y != 1) {
			return 1;
		}
	}
	{
		using a = q32::val<i32::val<1>, i32::val<1>>;
		using b = q32::val<i32::val<2>, i32::val<1>>;

		using c = q32::div_t<a, b>;

		auto x = c::x::v;
		auto y = c::y::v;
		if (x != 1 || y != 2) {
			return 1;
		}
	}

	return 0;
}

int test_simplify_q32() {
	using A = q32::val<i32::val<2>, i32::val<2>>;
	using B = q32::val<i32::val<1>, i32::val<1>>;
	using C = q32::val<i32::val<-1>, i32::val<-1>>;
	using D = q32::val<i32::val<1>, i32::val<-2>>;
	if (!q32::eq_v<q32::simplify_t<A>, q32::one>) {
		return 1;
	}
	if (!q32::eq_v<q32::simplify_t<B>, q32::one>) {
		return 1;
	}
	if (!q32::eq_v<q32::simplify_t<C>, q32::one>) {
		return 1;
	}
	if (!q32::eq_v<q32::simplify_t<D>, q32::val<i32::val<-1>, i32::val<2>>>) {
		return 1;
	}
	if(!q32::eq_v<q32::simplify_t<q32::sub_t<q32::zero, q32::zero>>, q32::zero>) {
		return 1;
	}

	return 0;
}

int test_eq_q32() {
	using A = q32::val<i32::val<2>, i32::val<2>>;
	using B = q32::val<i32::val<1>, i32::val<1>>;
	using C = q32::val<i32::val<-1>, i32::val<-1>>;
	if (!q32::eq_v<A, B>) {
		return 1;
	}
	if (!q32::eq_v<A, C>) {
		return 1;
	}
	return 0;
}

int test_fraction_field_of_fraction_field () {
	using qq32 = FractionField<q32>;
	if (!std::is_same<q32, qq32>::value) {
		return 1;
	}

	return 0;
}

int test_quotient_ring_is_z2z() {
	using QQ = Quotient<i32, i32::inject_constant_t<2>>;
	using two = QQ::inject_constant_t<2>;
	using three = QQ::inject_constant_t<3>;
	using one = QQ::inject_constant_t<1>;
	using zero = QQ::zero;
	if (!QQ::eq_v<two, zero>) {
		return 1;
	}
	if (!QQ::eq_v<one, three>) {
		return 1;
	}

	return 0;
}

int test_instanciate_F4() {
	using F2 = zpz<2>;
	using PF2 = polynomial<F2>;
	using F4 = Quotient<PF2, ConwayPolynomial<2, 2>::type>;
	using zero = F4::zero;
	using one = F4::one;
	using phi = F4::inject_ring_t<PF2::X>;
	using phi2 = F4::inject_ring_t<PF2::template val<F2::inject_constant_t<1>, F2::inject_constant_t<1>>>;
	// check that elements are different
	if (F4::eq_v<one, zero> || F4::eq_v<phi, zero> || F4::eq_v<phi2, zero> ||
		F4::eq_v<phi, one> || F4::eq_v<phi2, one> ||
		F4::eq_v<phi, phi2>) {
		return 1;
	}
	// one + phi = phi2
	if (!F4::eq_v<phi2, F4::template add_t<one, phi>>) {
		return 1;
	}
	// one + phi2 = phi
	if (!F4::eq_v<phi, F4::template add_t<one, phi2>>) {
		return 1;
	}
	// phi * phi = phi2
	if (!F4::eq_v<phi2, F4::template mul_t<phi, phi>>) {
		return 1;
	}
	// phi2 * phi2 = phi
	if (!F4::eq_v<phi, F4::template mul_t<phi2, phi2>>) {
		return 1;
	}
	return 0;
}

int test_instanciate_large_finite_field() {
	using F17 = zpz<17>;
	using PF17 = polynomial<F17>;
	using F_17_8 = Quotient<PF17, ConwayPolynomial<17, 8>::type>;
	using x = F_17_8::inject_ring_t<PF17::X>;
	return 0;
}

int test_factorial() {
	constexpr float x = factorial_v<i32, 3>;
	if (x != 6.0f) {
		return 1;
	}

	constexpr size_t y = factorial_v<i32, 0>;
	if (y != 1) {
		return 1;
	}

	return 0;
}

int test_combination() {
	constexpr int x = combination_v<i32, 2, 4>;
	if (x != 6) {
		return 1;
	}
	constexpr int y = combination_v<i32, 0, 4>;
	if (y != 1) {
		return 1;
	}
	constexpr int z = combination_v<i32, 1, 4>;
	if (z != 4) {
		return 1;
	}

	constexpr int zz = combination_v<i32, 3, 4>;
	if (zz != 4) {
		return 1;
	}
	return 0;
}

int test_bernouilli() {
	constexpr float b0 = bernouilli_v<float, i32, 0>;
	if (b0 != 1.0f) {
		return 1;
	}
	constexpr float b1 = bernouilli_v<float, i32, 1>;
	if (b1 != -0.5f) {
		return 1;
	}
	using B2 = bernouilli_t<i32, 2>;
	if (B2::x::v != 1) {
		return 1;
	}
	if (B2::y::v != 6) {
		return 1;
	}
	constexpr double b3 = bernouilli_v<double, i32, 3>;
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
	if (!is_prime<7883>::value) {
		return 1;
	}
	if (is_prime<7884>::value) {
		return 1;
	}
	if (!is_prime<7919>::value) {
		return 1;
	}
	if (is_prime<7920>::value) {
		return 1;
	}
	if (!is_prime<7927>::value) {
		return 1;
	}
	if (!is_prime<1000423>::value) {
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
	constexpr int a0 = internal::alternate<i32, 0>::value;
	if (a0 != 1) {
		return 1;
	}
	constexpr int a1 = internal::alternate<i32, 1>::value;
	if (a1 != -1) {
		return 1;
	}
	constexpr int a2 = internal::alternate<i32, 2>::value;
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
	static_assert(aerobus::IsRing<aerobus::q32>);
	static_assert(aerobus::IsField<aerobus::q32>);
	static_assert(aerobus::IsRing<aerobus::polynomial<i32>>);
	static_assert(aerobus::IsField<aerobus::FractionField<aerobus::polynomial<i32>>>);
	static_assert(IsRing<bigint>);
	return 0;
}

int test_continued_fraction() {
	// A001203
	constexpr double A_PI = PI_fraction::val;
	if(::fabs(A_PI - M_PI) > 0) {
		return 1;
	}

	// A003417
	constexpr double A_E = E_fraction::val;
	if(::fabs(A_E - M_E) > 0) {
		return 1;
	}

	constexpr double A_SQRT2 = SQRT2_fraction::val;
	if(::fabs(A_SQRT2 - M_SQRT2) > 0) {
		return 1;
	}

	constexpr double A_SQRT3 = SQRT3_fraction::val;
	if(::fabs(A_SQRT3 - 1.7320508075688772935) > 0) {
		::printf("%.17g", A_SQRT3);
		return 1;
	}

	return 0;
}

int test_chebyshev() {
	using T4 = aerobus::chebyshev_T<4>;

	if(T4::degree != 4) {
		return 1;
	}
	if (T4::template coeff_at_t<4>::v != 8) {
		return 0;
	}
	
	if (T4::template coeff_at_t<3>::v != 0) {
		return 0;
	}
	if (T4::template coeff_at_t<2>::v != -8) {
		return 0;
	}
	if (T4::template coeff_at_t<1>::v != 0) {
		return 0;
	}
	if (T4::template coeff_at_t<0>::v != 1) {
		return 0;
	}

	using U4 = chebyshev_U<4>;
	if(U4::degree != 4) {
		return 1;
	}
	if (U4::template coeff_at_t<4>::v != 16) {
		return 0;
	}
	if (U4::template coeff_at_t<3>::v != 0) {
		return 0;
	}
	if (U4::template coeff_at_t<2>::v != -12) {
		return 0;
	}
	if (U4::template coeff_at_t<1>::v != 0) {
		return 0;
	}
	if (U4::template coeff_at_t<0>::v != 1) {
		return 0;
	}


	return 0;
}

int test_value_at() {
	{
		constexpr uint16_t x = aerobus::internal::value_at<0, 1, 2, 3>::value;
		if(x != 1) {
			return 1;
		}
	}
	{
		constexpr uint16_t x = aerobus::internal::value_at<1, 1, 2, 3>::value;
		if(x != 2) {
			return 1;
		}
	}
	{
		constexpr uint16_t x = aerobus::internal::value_at<2, 1, 2, 3>::value;
		if(x != 3) {
			return 1;
		}
	}
	{
		constexpr uint16_t x = aerobus::internal::value_at<0, 1, 2, 3>::value;
		if(x != 1) {
			return 1;
		}
	}
	return 0;
}

int test_bigint_digit_at() {
	using x = bigint_pos<1, 2>;
	constexpr uint16_t x0 = x::template digit_at<0>::value;
	if(x0 != 2) {
		return 1;
	}
	constexpr uint16_t x1 = x::template digit_at<1>::value;
	if(x1 != 1) {
		return 1;
	}
	constexpr uint16_t x2 = x::template digit_at<2>::value;
	if(x2 != 0) {
		return 1;
	}

	return 0;
}

int test_bigint_minus() {
	using I1 = bigint_pos<1, 2>;
	using I2 = bigint_neg<1, 2>;
	using MZ = bigint_neg<0>;
	if(!bigint::eq_v<I1, typename I2::minus_t>) {
		return 1;
	}
	if(!bigint::eq_v<I2, typename I1::minus_t>) {
		return 1;
	}
	if(!bigint::eq_v<typename bigint::zero, MZ>) {
		return 1;
	}

	return 0;
}

int test_bigint_add() {
	using I1 = bigint_pos<UINT32_MAX, UINT32_MAX>;
	using I2 = bigint::one;
	using S = bigint::add_t<I1, I2>;
	using S2 = bigint::add_t<I1, typename bigint::zero>;
	using S3 = bigint::add_t<bigint::zero, I2>;
	if(!bigint::eq_v<S2, I1> || !bigint::eq_v<I2, S3>) {
		return 1;
	}
	if(S::digits != 3) {
		return 1;
	}
	if(S::digit_at<0>::value != 0) {
		return 0;
	}
	if(S::digit_at<1>::value != 0) {
		return 0;
	}
	if(S::digit_at<2>::value != 1) {
		return 0;
	}
	// 22 + -11 = 11
	{
		using X = bigint_pos<2, 2>;
		using Y = bigint_neg<1, 1>;
		using S4 = bigint::add_t<X, Y>;
		if(!bigint::eq_v<S4, typename Y::minus_t>) {
			printf("%s\n", S4::to_string().c_str());
			return 1;
		}
	}

	// -22 + -11 = -33
	{
		using X = bigint_neg<2, 2>;
		using Y = bigint_neg<1, 1>;
		using expected = bigint_neg<3, 3>;
		using S4 = bigint::add_t<X, Y>;
		if(!bigint::eq_v<S4, expected>) {
			return 1;
		}
	}

	// -22 + 11 = -11
	{
		using X = bigint_neg<2, 2>;
		using Y = bigint_pos<1, 1>;
		using expected = bigint_neg<1, 1>;
		using S4 = bigint::add_t<X, Y>;
		if(!bigint::eq_v<S4, expected>) {
			printf("%s\n", S4::to_string().c_str());
			return 1;
		}
	}

	return 0;
}

int test_bigint_sub() {
	// 100 - 99 = 1
	{
		using b = bigint_pos<UINT32_MAX, UINT32_MAX>;
		using c = bigint_pos<1, 0, 0>;
		using d = bigint::sub_t<c, b>;
		if(!bigint::eq_v<bigint::one, d>) {
			return 1;
		}
	}
	// 100 - 1 = 99
	{
		using b = bigint_pos<UINT32_MAX, UINT32_MAX>;
		using c = bigint_pos<1, 0, 0>;
		using d = bigint::sub_t<c, bigint::one>;
		if(!bigint::eq_v<b, d>) {
			return 1;
		}
	}
	// 22 - 11 = 11
	{
		using a = bigint_pos<2, 2>;
		using b = bigint_pos<1, 1>;
		using c = bigint::sub_t<a, b>;
		if(!bigint::eq_v<b, c>) {
			return 1;
		}
	}
	// 1 - 2 = -1
	{
		using a = bigint_pos<1>;
		using b = bigint_pos<2>;
		using expected = bigint_neg<1>;
		using c = bigint::sub_t<a, b>;
		if(!bigint::eq_v<expected, c>) {
			return 1;
		}
	}
	// 2 - -1 = 3
	{
		using a = bigint_pos<2>;
		using b = bigint_neg<1>;
		using expected = bigint_pos<3>;
		using c = bigint::sub_t<a, b>;
		if(!bigint::eq_v<expected, c>) {
			return 1;
		}
	}
	// -2 - -1 = -1
	{
		using a = bigint_neg<2>;
		using b = bigint_neg<1>;
		using expected = bigint_neg<1>;
		using c = bigint::sub_t<a, b>;
		if(!bigint::eq_v<expected, c>) {
			return 1;
		}
	}
	// 99 - -1 = -1
	{
		using a = bigint_pos<UINT32_MAX, UINT32_MAX>;
		using b = bigint_neg<1>;
		using expected = bigint_pos<1, 0, 0>;
		using c = bigint::sub_t<a, b>;
		if(!bigint::eq_v<expected, c>) {
			return 1;
		}
	}
	// 99 - 0 = 00
	{
		using a = bigint_pos<UINT32_MAX, UINT32_MAX>;
		using c = bigint::sub_t<a, bigint::zero>;
		if(!bigint::eq_v<a, c>) {
			return 1;
		}
	}
	// 12 - 11 = 1 (must be simplified)
	{
		using a = bigint_pos<1, 2>;
		using b = bigint_pos<1, 1>;
		using c = bigint::sub_t<a, b>;
		if(!bigint::eq_v<bigint::one, c>) {
			return 1;
		}
	}

	return 0; 
}

int test_bigint_mul() {
	{
		using A = bigint_pos<1, 2>;
		using B = bigint::mul_t<A, bigint::one>;
		if(!bigint::eq_v<A, B>) {
			return 1;
		}
	}
	{
		using A = bigint_pos<1, 2>;
		using B = bigint::mul_t<A, bigint::zero>;
		if(!bigint::eq_v<bigint::zero, B>) {
			return 1;
		}
	}
	{
		using A = bigint_pos<1, 2>;
		using B = bigint_pos<2>;
		using PROD = bigint::mul_t<A, B>;
		using expected = bigint_pos<2, 4>;
		if(!bigint::eq_v<expected, PROD>) {
			printf("%s\n", PROD::to_string().c_str());
			return 1;
		}
	}
	// 9 * 9 = 81
	{
		using A = bigint_pos<UINT32_MAX>;
		using B = bigint_pos<UINT32_MAX>;
		using PROD = bigint::mul_t<A, B>;
		using expected = bigint_pos<UINT32_MAX - 1, 1>;
		if(!bigint::eq_v<expected, PROD>) {
			printf("%s\n", PROD::to_string().c_str());
			return 1;
		}
	}

	// 9 * -9 = -81
	{
		using A = bigint_pos<UINT32_MAX>;
		using B = bigint_neg<UINT32_MAX>;
		using PROD = bigint::mul_t<A, B>;
		using expected = bigint_neg<UINT32_MAX - 1, 1>;
		if(!bigint::eq_v<expected, PROD>) {
			printf("%s\n", PROD::to_string().c_str());
			return 1;
		}
	}	

	// -9 * -9 = 81
	{
		using A = bigint_neg<UINT32_MAX>;
		using B = bigint_neg<UINT32_MAX>;
		using PROD = bigint::mul_t<A, B>;
		using expected = bigint_pos<UINT32_MAX - 1, 1>;
		if(!bigint::eq_v<expected, PROD>) {
			printf("%s\n", PROD::to_string().c_str());
			return 1;
		}
	}	
	// 99 * 9 = 891
	{
		using A = bigint_pos<UINT32_MAX, UINT32_MAX>;
		using B = bigint_pos<UINT32_MAX>;
		using PROD = bigint::mul_t<A, B>;
		using expected = bigint_pos<UINT32_MAX - 1, UINT32_MAX, 1>;
		if(!bigint::eq_v<expected, PROD>) {
			printf("%s\n", PROD::to_string().c_str());
			return 1;
		}
	}
	// 999 * 99 = 98901
	{
		using A = bigint_pos<UINT32_MAX, UINT32_MAX, UINT32_MAX>;
		using B = bigint_pos<UINT32_MAX, UINT32_MAX>;
		using PROD = bigint::mul_t<A, B>;
		using expected = bigint_pos<UINT32_MAX, UINT32_MAX - 1, UINT32_MAX, 0, 1>;
		if(!bigint::eq_v<expected, PROD>) {
			printf("%s\n", PROD::to_string().c_str());
			return 1;
		}
	}

	return 0;
}

int test_bigint_eq() {
	using a = bigint::val<bigint::signs::positive, 1, 2, 3>;
	using b = bigint::val<bigint::signs::negative, 1, 2, 3>;
	using c = bigint::val<bigint::signs::positive, 1, 2, 1>;
	if(bigint::eq_v<a, b>) {
		return 1;
	}
	if(!bigint::eq_v<a, a>) {
		return 1;
	}
	if(bigint::eq_v<a, c>) {
		return 1;
	}

	return 0;
}

int test_bigint_pos() {
	using a = bigint::zero;
	using b = bigint::val<bigint::signs::positive, 1>;
	using c = bigint::val<bigint::signs::negative, 1>;
	using d = bigint::val<bigint::signs::negative, 0>;

	if(bigint::pos_v<a>) {
		return 1;
	}
	if(!bigint::pos_v<b>) {
		return 1;
	}
	if(bigint::pos_v<c>) {
		return 1;
	}
	if(bigint::pos_v<d>) {
		return 1;
	}

	return 0;
}

int test_bigint_gt() {
	using a = bigint::zero;
	using aa = bigint_neg<0>; // negative zero -- should not exist
	using b = bigint_pos<1, 2, 3>;
	using c = bigint_neg<1>;
	using d = bigint_neg<2>;
	using e = bigint_neg<1, 2>;
	using f = bigint::one;
	// zero
	if(bigint::gt_v<a, b> || bigint::gt_v<a, aa> || bigint::gt_v<a, aa> || bigint::gt_v<aa, a>) {
		return 1;
	}
	// zero
	if(!bigint::gt_v<a, c> || !bigint::gt_v<a, e>) {
		return 1;
	}
	// normal positive numbers
	if(!bigint::gt_v<b, a> || !bigint::gt_v<b, f>) {
		return 1;
	}
	// normal negative numbers
	if(!bigint::gt_v<c, d> || !bigint::gt_v<c, e> || !bigint::gt_v<d, e>) {
		return 1;
	}

	return 0;
}

int test_bigint_shift_left() {
	{
		using X = bigint_pos<2>;
		using Y = bigint::shift_left_t<X, 2>;
		using expected = bigint_pos<2, 0, 0>;
		if(!bigint::eq_v<Y, expected>) {
			return 1;
		}
	}
	{
		using X = bigint_pos<2>;
		using Y = bigint::shift_left_t<X, 0>;
		if(!bigint::eq_v<X, Y>) {
			return 1;
		}
	}
	{
		using X = bigint_pos<2, 3>;
		using Y = bigint::shift_left_t<X, 2>;
		using expected = bigint_pos<2, 3, 0, 0>;
		if(!bigint::eq_v<Y, expected>) {
			return 1;
		}
	}

	return 0;
}

int test_bigint_shift_right() {
	{
		using A = bigint_pos<1, 2, 3>;
		using B = bigint::shift_right_t<A, 1>;
		if (!bigint::eq_v < B, bigint_pos<1, 2>>) {
			printf("%s\n", B::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<1, 2, 3>;
		using B = bigint::shift_right_t<A, 2>;
		if (!bigint::eq_v < B, bigint_pos<1>>) {
			printf("%s\n", B::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<1, 2, 3>;
		using B = bigint::shift_right_t<A, 3>;
		if (!bigint::eq_v<B, typename bigint::zero>) {
			printf("%s\n", B::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<1, 2, 3>;
		using B = bigint::shift_right_t<A, 0>;
		if (!bigint::eq_v<B, A>) {
			printf("%s\n", B::to_string().c_str());
			return 1;
		}
	}
	return 0;
}

int test_bigint_div2() {
	{
		using A = bigint_pos<12>;
		using B = bigint::div_2_t<A>;
		if (!bigint::eq_v<B, bigint_pos<6>>) {
			return 1;
		}
	} 
	{
		using A = bigint_pos<1, 2>;
		using B = bigint::div_2_t<A>;
		using expected = bigint_pos<0x8000'0001U>;
		if (!bigint::eq_v<B, expected>) {
			printf("%s\n", B::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<UINT32_MAX, 2>;
		using B = bigint::div_2_t<A>;
		using C = bigint::mul_t<B, bigint_pos<2>>;
		if (!bigint::eq_v<A, C>) {
			printf("%s\n", B::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<UINT32_MAX, UINT32_MAX, UINT32_MAX-1>;
		using B = bigint::div_2_t<A>;
		using C = bigint::mul_t<B, bigint_pos<2>>;
		if (!bigint::eq_v<A, C>) {
			printf("%s\n", B::to_string().c_str());
			return 1;
		}
	}

	return 0;
}

int test_bigint_floor() {
	{
		using A = bigint_pos<12>;
		using B = bigint_pos<5>;
		using C = bigint::floor_t<A, B>;
		if (!bigint::eq_v<C, bigint_pos<2>>) {
			return 1;
		}
	}
	{
		using A = bigint_pos<UINT32_MAX>;
		using B = bigint_pos<5>;
		using C = bigint::floor_t<A, B>;
		using C1 = bigint::add_t<C, bigint::one>;
		using mul = bigint::mul_t<C, B>;
		using mul1 = bigint::mul_t<C1, B>;

		if (!bigint::ge_v<A, mul> || !bigint::gt_v<mul1, A>) {
			printf("%s\n", C::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<UINT32_MAX, 1>;
		using B = bigint_pos<2>;
		using C = bigint::floor_t<A, B>;
		using C1 = bigint::add_t<C, bigint::one>;
		using mul = bigint::mul_t<C, B>;
		using mul1 = bigint::mul_t<C1, B>;

		if (!bigint::ge_v<A, mul> || !bigint::gt_v<mul1, A>) {
			printf("floor = %s\n", C::to_string().c_str());
			printf("original = %s\n", A::to_string().c_str());
			printf("mul = %s\n", mul::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<UINT32_MAX, 1, 7, 0xFF23140A>;
		using B = bigint_pos<UINT32_MAX, 1>;
		using C = bigint::floor_t<A, B>;
		using C1 = bigint::add_t<C, bigint::one>;
		using mul = bigint::mul_t<C, B>;
		using mul1 = bigint::mul_t<C1, B>;

		if (!bigint::ge_v<A, mul> || !bigint::gt_v<mul1, A>) {
			printf("floor = %s\n", C::to_string().c_str());
			printf("original = %s\n", A::to_string().c_str());
			printf("mul = %s\n", mul::to_string().c_str());
			return 1;
		}
	}

	{
		using A = bigint_pos<128, 41>;
		using B = bigint_pos<25, 1857>;
		using C = bigint::floor_t<A, B>;
		using C1 = bigint::add_t<C, bigint::one>;
		using mul = bigint::mul_t<C, B>;
		using mul1 = bigint::mul_t<C1, B>;

		if (!bigint::ge_v<A, mul> || !bigint::gt_v<mul1, A>) {
			printf("floor = %s\n", C::to_string().c_str());
			printf("original = %s\n", A::to_string().c_str());
			printf("mul = %s\n", mul::to_string().c_str());
			return 1;
		}
	}

	return 0;
}

int test_bigint_div() {
	{
		using A = bigint_pos<12>;
		using B = bigint_pos<6>;
		using DIV = bigint::div_t<A, B>;
		using REM = bigint::mod_t<A, B>;
		using QEXP = bigint_pos<2>;
		using REXP = bigint_pos<0>;
		if(!bigint::eq_v<QEXP, DIV> || !bigint::eq_v<REXP, REM>) {
			printf("Q = %s\nR = %s\n", DIV::to_string().c_str(), REM::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<12, 12>;
		using B = bigint_pos<6>;
		using DIV = bigint::div_t<A, B>;
		using REM = bigint::mod_t<A, B>;
		using QEXP = bigint_pos<2, 2>;
		using REXP = bigint_pos<0>;
		if(!bigint::eq_v<QEXP, DIV> || !bigint::eq_v<REXP, REM>) {
			printf("Q = %s\nR = %s\n", DIV::to_string().c_str(), REM::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<12, 13>;
		using B = bigint_pos<6>;
		using DIV = bigint::div_t<A, B>;
		using REM = bigint::mod_t<A, B>;
		using QEXP = bigint_pos<2, 2>;
		using REXP = bigint_pos<1>;
		if(!bigint::eq_v<QEXP, DIV> || !bigint::eq_v<REXP, REM>) {
			printf("Q = %s\nR = %s\n", DIV::to_string().c_str(), REM::to_string().c_str());
			printf("QEXP = %s\nR = %s\n", QEXP::to_string().c_str(), REXP::to_string().c_str());
			return 1;
		}
	}

	{
		using A = bigint_pos<12, 13>;
		using B = bigint_pos<12, 13>;
		using DIV = bigint::div_t<A, B>;
		using REM = bigint::mod_t<A, B>;
		using QEXP = bigint_pos<1>;
		using REXP = bigint_pos<0>;
		if(!bigint::eq_v<QEXP, DIV> || !bigint::eq_v<REXP, REM>) {
			printf("Q = %s\nR = %s\n", DIV::to_string().c_str(), REM::to_string().c_str());
			printf("QEXP = %s\nR = %s\n", QEXP::to_string().c_str(), REXP::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<12, 13>;
		using B = bigint_pos<2>;
		using DIV = bigint::div_t<A, B>;
		using REM = bigint::mod_t<A, B>;
		using QEXP = bigint_pos<6, 6>;
		using REXP = bigint_pos<1>;
		if(!bigint::eq_v<QEXP, DIV> || !bigint::eq_v<REXP, REM>) {
			printf("Q = %s\nR = %s\n", DIV::to_string().c_str(), REM::to_string().c_str());
			printf("QEXP = %s\nR = %s\n", QEXP::to_string().c_str(), REXP::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<12, 12>;
		using B = bigint_pos<9>;
		using DIV = bigint::div_t<A, B>;
		using REM = bigint::mod_t<A, B>;

		using AA = bigint::fma_t<B, DIV, REM>;
		if(!bigint::eq_v<AA, A>) {
			printf("Q = %s\nR = %s\n", DIV::to_string().c_str(), REM::to_string().c_str());
			return 1;
		}
	}
	{
		using A = bigint_pos<128, 41, 17899>;
		using B = bigint_pos<25, 1857>;
		using D = bigint::div_t<A, B>;
		using R = bigint::mod_t<A, B>;
		using E = bigint::fma_t<B, D, R>;
		if(!bigint::eq_v<E, A>) {
			printf("D = %s\nR = %s\n", D::to_string().c_str(), R::to_string().c_str());
			return 1;
		}
	}

	return 0;
}

int test_bigint_gcd() {
	{
		using A = bigint_pos<12>;
		using B = bigint_pos<3>;
		using G = bigint::gcd_t<A, B>;
		if(!bigint::eq_v<bigint_pos<3>, G>) {
			return 1;
		}
	}
	{
		using A = bigint_pos<12, 12>;
		using B = bigint_pos<9>;
		using G = bigint::gcd_t<A, B>;
		if(!bigint::eq_v<bigint_pos<3>, G>) {
			return 1;
		}
	}

	return 0;
}

static uint32_t ok_count = 0;
static uint32_t fail_count = 0;
#define RUN_TEST(test_name) \
	if (test_name() != 0) { \
		printf("%s failed\n", #test_name); \
		fail_count += 1; \
	} else { \
		ok_count += 1; \
	}

int main(int argc, char* argv[]) {
	RUN_TEST(test_chebyshev);
	RUN_TEST(test_continued_fraction)
	RUN_TEST(test_concept_ring)
	RUN_TEST(test_type_list)
	RUN_TEST(test_type_at) 
	RUN_TEST(test_poly_simplify)
	RUN_TEST(test_coeff_at)
	RUN_TEST(test_poly_add)
	RUN_TEST(test_poly_sub)
	RUN_TEST(test_poly_derive)
	RUN_TEST(test_poly_eq)
	RUN_TEST(test_gcd)
	RUN_TEST(test_poly_mul)
	RUN_TEST(test_poly_to_string)
	RUN_TEST(test_monomial)
	RUN_TEST(test_poly_gcd)
	RUN_TEST(test_poly_eval)
	RUN_TEST(test_add_q32)
	RUN_TEST(test_gt_q32)
	RUN_TEST(test_is_zero_q32)
	RUN_TEST(test_fraction_field_eval)
	RUN_TEST(test_sub_q32)
	RUN_TEST(test_mul_q32)
	RUN_TEST(test_div_q32)
	RUN_TEST(test_eq_q32)
	RUN_TEST(test_fraction_field_of_fraction_field)
	RUN_TEST(test_quotient_ring_is_z2z)
	RUN_TEST(test_instanciate_F4)
	RUN_TEST(test_instanciate_large_finite_field)
	RUN_TEST(test_factorial)
	RUN_TEST(test_combination)
	RUN_TEST(test_bernouilli)
	RUN_TEST(test_alternate)
	RUN_TEST(test_exp)
	RUN_TEST(test_is_prime)
	RUN_TEST(test_zpz)
	RUN_TEST(test_value_at)
	RUN_TEST(test_bigint_digit_at)
	RUN_TEST(test_bigint_minus)
	RUN_TEST(test_bigint_add)
	RUN_TEST(test_bigint_sub)
	RUN_TEST(test_bigint_eq)
	RUN_TEST(test_bigint_pos)
	RUN_TEST(test_bigint_gt)
	RUN_TEST(test_bigint_shift_left)
	RUN_TEST(test_bigint_shift_right)
	RUN_TEST(test_bigint_mul)
	RUN_TEST(test_bigint_div2)
	RUN_TEST(test_bigint_floor)
	RUN_TEST(test_bigint_div)
	RUN_TEST(test_bigint_gcd);

	printf("%d/%d tests passed\n", ok_count, ok_count + fail_count);
	if(fail_count > 0) {
		return 1;
	}

	if(fail_count == 0) {
		printf("All tests passed\n");
	}
	return 0;
}