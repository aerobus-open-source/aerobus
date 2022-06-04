#include <cstdio>

#include "lib.h"

int test_type_at() {
	if (!std::is_same<type_at_t<0, float, int, long>, float>::value) {
		return 1;
	}
	if (!std::is_same<type_at_t<1, float, int, long>, int>::value) {
		return 1;
	}
	if (!std::is_same<type_at_t<2, float, int, long>, long>::value) {
		return 1;
	}

	return 0;
}

int test_poly_simplify() {
	using poly1 = polynomial<i32>::val<i32::val<0>, i32::val<1>, i32::val<2>>;
	using simplified1 = poly_simplify_t<i32, poly1>;
	using expected1 = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
	if (!std::is_same<expected1, simplified1>::value) {
		return 1;
	}

	using poly2 = polynomial<i32>::val<i32::val<12>>;
	using simplified2 = poly_simplify_t<i32, poly2>;
	if (!std::is_same<poly2, simplified2>::value) {
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

int test_poly_add_at() {
	{
		using P1 = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
		using P2 = polynomial<i32>::val<i32::val<2>, i32::val<3>>;
		using add0 = poly_add_at_t<i32, P1, P2, 0>;
		using add1 = poly_add_at_t<i32, P1, P2, 1>;
		using add2 = poly_add_at_t<i32, P1, P2, 2>;
		{
			auto expected = 5;
			auto actual = add0::v;
			if (expected != actual) {
				printf("expected %d -- got %d\n", expected, actual);
				return 1;
			}
		} 
		{
			auto expected = 3;
			auto actual = add1::v;
			if (expected != actual) {
				printf("expected %d -- got %d\n", expected, actual);
				return 1;
			}
		}
		{
			auto expected = 0;
			auto actual = add2::v;
			if (expected != actual) {
				printf("expected %d -- got %d\n", expected, actual);
				return 1;
			}
		}
	}

	return 0;
}


template<typename... coeffs>
using IX = polynomial<i32>::val<coeffs...>;
template<int32_t x>
using Int = i32::val<x>;

template<typename P1, typename P2>
using add_ix = poly_add_t<i32, P1, P2>;
template<typename P1, typename P2>
using sub_ix = poly_sub_t<i32, P1, P2>;

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
		using A = poly_add_t<i32, P1, P2>;
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
		using A = poly_sub_t<i32, P2, P1>;
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
	using A = i32::gcd_t<i32::val<6>, i32::val<12>>;
	if (A::v != 6) {
		return 1;
	}
	using B = i32::gcd_t<i32::val<12>, i32::val<6>>;
	if (B::v != 6) {
		return 1;
	}

	using C = i32::gcd_t<i32::val<3>, i32::val<5>>;
	if (C::v != 1) {
		return 1;
	}

	return 0;
}

int main(int argc, char* argv[]) {
	if (test_type_at() != 0) {
		return 1;
	}
	if (test_poly_simplify() != 0) {
		return 1;
	}
	if (test_coeff_at() != 0) {
		return 1;
	}
	if (test_poly_add_at() != 0) {
		return 1;
	}

	if (test_poly_add() != 0) {
		return 1;
	}
	if (test_poly_sub() != 0) {
		return 1;
	}
	if (test_poly_eq() != 0) {
		return 1;
	}
	if (test_gcd() != 0) {
		return 1;
	}

	return 0;
}