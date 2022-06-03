#include <cstdio>
#include <typeinfo>

#include "lib.h"

int test_type_at() {
	const char* a = typeid(type_at_t<0, float, int, long>).name();
	const char* b = typeid(type_at_t<1, float, int, long>).name();
	const char* c = typeid(type_at_t<2, float, int, long>).name();
	if (strcmp(a, "float") != 0) {
		printf("expected float, got %s\n", a);
		return 1;
	}
	if (strcmp(b, "int") != 0) {
		printf("expected int, got %s\n", a);
		return 1;
	}
	if (strcmp(c, "long") != 0) {
		printf("expected long, got %s\n", a);
		return 1;
	}

	return 0;
}

int test_poly_simplify() {
	using poly1 = polynomial<i32>::val<i32::val<0>, i32::val<1>, i32::val<2>>;
	using simplified1 = poly_simplify_t<i32, poly1>;
	if (strcmp(typeid(simplified1).name(), "struct polynomial<struct i32>::val<struct i32::val<1>,struct i32::val<2> >") != 0) {
		printf("got %s\n", typeid(simplified1).name());
		return 1;
	}

	using poly2 = polynomial<i32>::val<i32::val<12>>;
	using simplified2 = poly_simplify_t<i32, poly2>;
	if (strcmp(typeid(simplified2).name(), "struct polynomial<struct i32>::val<struct i32::val<12> >") != 0) {
		printf("got %s\n", typeid(simplified2).name());
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

	if (strcmp(typeid(at0).name(), typeid(i32::val<1>).name()) != 0) {
		return 1;
	}
	if (strcmp(typeid(at1).name(), typeid(i32::val<2>).name()) != 0) {
		return 1;
	}
	if (strcmp(typeid(at2).name(), typeid(i32::val<3>).name()) != 0) {
		return 1;
	}
	if (strcmp(typeid(at3).name(), typeid(i32::zero).name()) != 0) {
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
		using A = add_ix<P1, P2>;

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
	// 2 + x
	/*using poly = polynomial<i32>::val<i32::val<0>, i32::val<1>, i32::val<2>>;
	using simplified = poly_simplify_t<i32, poly>;
	auto x = typeid(simplified).name();*/

	//using at = simplified::coeff_at<1>;
	//auto y = typeid(at).name();
	//printf("%s\n", y); // struct i32::val<1>

	//printf("%s\n", x); // struct polynomial<struct i32>::val<struct i32::val<1>,struct i32::val<1> >
}