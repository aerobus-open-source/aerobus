#include <cstdint>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <algorithm>

template <int64_t i, typename T, typename... Ts>
struct type_at
{
	static_assert(i < sizeof...(Ts) + 1, "index out of range");
	using type = typename type_at<i - 1, Ts...>::type;
};

template <typename T, typename... Ts> struct type_at<0, T, Ts...> {
	using type = T;
};

template <size_t i, typename... Ts>
using type_at_t = typename type_at<i, Ts...>::type;

template <std::size_t ... Is>
constexpr auto index_sequence_reverse(std::index_sequence<Is...> const&)
-> decltype(std::index_sequence<sizeof...(Is) - 1U - Is...>{});

template <std::size_t N>
using make_index_sequence_reverse
= decltype(index_sequence_reverse(std::make_index_sequence<N>{}));

template<typename Ring>
struct gcd {
	template<typename a, typename b, typename E = void>
	struct gcd_helper {};

	template<typename a, typename b>
	struct gcd_helper<a, b, std::enable_if_t<
		std::is_same<
		typename Ring::template eq_t<a, b>,
		std::true_type
		>::value
	>
	> {
		using type = a;
	};

	template<typename a, typename b>
	struct gcd_helper<a, b, std::enable_if_t<
		std::is_same<a, typename Ring::zero>::value || std::is_same<b, typename Ring::zero>::value>> {
		using type = Ring::zero;
	};

	template<typename a, typename b>
	struct gcd_helper<a, b, std::enable_if_t<
		std::is_same<
		typename Ring::template gt_t<a, b>,
		std::true_type
		>::value &&
		!std::is_same<a, typename Ring::zero>::value &&
		!std::is_same<b, typename Ring::zero>::value>
	> {
		using type = typename gcd_helper<typename Ring::template sub_t<a, b>, b>::type;
	};

	template<typename a, typename b>
	struct gcd_helper<a, b, std::enable_if_t<
		std::is_same<
		typename Ring::template lt_t<a, b>,
		std::true_type
		>::value &&
		!std::is_same<a, typename Ring::zero>::value &&
		!std::is_same<b, typename Ring::zero>::value>
	> {
		using type = typename gcd_helper<a, typename Ring::template sub_t<b, a>>::type;
	};

	template<typename a, typename b>
	using type = typename gcd<Ring>::template gcd_helper<a, b>::type;
};

struct i32 {
	template<int32_t x>
	struct val {
		static constexpr int32_t v = x;
	};

	using zero = val<0>;
	using one = val<1>;
	static constexpr bool is_field = false;

private:
	template<typename v1, typename v2>
	struct add {
		using type = val<v1::v + v2::v>;
	};

	template<typename v1, typename v2>
	struct sub {
		using type = val<v1::v - v2::v>;
	};

	template<typename v1, typename v2>
	struct mul {
		using type = val<v1::v* v2::v>;
	};

	template<typename v1, typename v2>
	struct div {
		using type = val<v1::v / v2::v>;
	};

	template<typename v1, typename v2>
	struct gt {
		using type = typename std::conditional<(v1::v > v2::v), std::true_type, std::false_type>::type;
	};

	template<typename v1, typename v2>
	struct lt {
		using type = typename std::conditional<(v1::v < v2::v), std::true_type, std::false_type>::type;
	};

	template<typename v1, typename v2>
	struct eq {
		using type = typename std::conditional<(v1::v == v2::v), std::true_type, std::false_type>::type;
	};

public:
	template<typename v1, typename v2>
	using add_t = typename add<v1, v2>::type;

	template<typename v1, typename v2>
	using sub_t = typename sub<v1, v2>::type;

	template<typename v1, typename v2>
	using mul_t = typename mul<v1, v2>::type;

	template<typename v1, typename v2>
	using div_t = typename div<v1, v2>::type;

	template<typename v1, typename v2>
	using gt_t = typename gt<v1, v2>::type;

	template<typename v1, typename v2>
	using lt_t = typename lt<v1, v2>::type;

	template<typename v1, typename v2>
	using eq_t = typename eq<v1, v2>::type;

	template<typename v1, typename v2>
	using gcd_t = typename gcd<i32>::template type<v1, v2>;
};

template<typename Ring, typename P, typename E = void>
struct poly_simplify;

template<typename Ring, typename P>
using poly_simplify_t = typename poly_simplify<Ring, P>::type;

template <typename Ring, typename P1, typename P2, typename I>
struct poly_add_low;

template <typename Ring, typename P1, typename P2, typename I>
struct poly_sub_low;

template <typename Ring, typename P1, typename P2, typename I>
struct poly_mul_low;

template<typename Ring, typename P1, typename P2>
using poly_add_t = poly_simplify_t<Ring, typename poly_add_low<
	Ring,
	P1,
	P2,
	make_index_sequence_reverse<
	std::max(P1::degree, P2::degree) + 1
	>>::type>;

template<typename Ring, typename P1, typename P2>
using poly_sub_t = poly_simplify_t<Ring, typename poly_sub_low<
	Ring,
	P1,
	P2,
	make_index_sequence_reverse<
	std::max(P1::degree, P2::degree) + 1
	>>::type>;

template<typename Ring, typename P1, typename P2>
using poly_mul_t = typename poly_mul_low<
	Ring,
	P1,
	P2,
	make_index_sequence_reverse<
	P1::degree + P2::degree + 1
	>>::type;

template<typename Ring, typename A, typename B, typename Q, typename R, typename E = void>
struct poly_div_helper;

// coeffN x^N + ...
template<typename Ring>
struct polynomial {
	static constexpr bool is_field = false;

	template<typename coeffN, typename... coeffs>
	struct val {
		static constexpr size_t degree = sizeof...(coeffs);
		using aN = coeffN;
		using strip = polynomial<Ring>::template val<coeffs...>;

		template<size_t index, typename E = void>
		struct coeff_at {};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index >= 0 && index <= sizeof...(coeffs))>> {
			using type = type_at_t<sizeof...(coeffs) - index, coeffN, coeffs...>;
		};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index < 0 || index > sizeof...(coeffs))>> {
			using type = typename Ring::zero;
		};

		template<size_t index>
		using coeff_at_t = typename coeff_at<index>::type;
	};

	// specialization for constants
	template<typename coeffN>
	struct val<coeffN> {
		static constexpr size_t degree = 0;
		using aN = coeffN;
		using strip = val<coeffN>;

		template<size_t index, typename E = void>
		struct coeff_at {};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index == 0)>> {
			using type = aN;
		};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index < 0 || index > 0)>> {
			using type = typename Ring::zero;
		};

		template<size_t index>
		using coeff_at_t = typename coeff_at<index>::type;
	};

private:
	template<typename v1, typename v2, typename E = void>
	struct eq_helper {};

	template<typename v1, typename v2>
	struct eq_helper<v1, v2, std::enable_if_t<v1::degree != v2::degree>> {
		using type = std::false_type;
	};

	template<typename v1, typename v2>
	struct eq_helper<v1, v2, std::enable_if_t<
		v1::degree == v2::degree &&
		v1::degree != 0 &&
		std::is_same<
		typename Ring::template eq_t<typename v1::aN, typename v2::aN>,
		std::true_type
		>::value
	>
	> {
		using type = typename eq_helper<typename v1::strip, typename v2::strip>::type;
	};

	template<typename v1, typename v2>
	struct eq_helper<v1, v2, std::enable_if_t<
		v1::degree == v1::degree &&
		v1::degree == 0
	>
	> {
		using type = typename Ring::template eq_t<typename v1::aN, typename v2::aN>;
	};

	template<typename v1, typename v2, typename E = void>
	struct lt_helper {};

	template<typename v1, typename v2>
	struct lt_helper<v1, v2, std::enable_if_t<(v1::degree < v2::degree)>> {
		using type = std::true_type;
	};

	template<typename v1, typename v2>
	struct lt_helper<v1, v2, std::enable_if_t<(v1::degree == v2::degree)>> {
		using type = typename Ring::template lt_t<typename v1::aN, typename v2::aN>;
	};

	template<typename v1, typename v2>
	struct lt_helper<v1, v2, std::enable_if_t<(v1::degree > v2::degree)>> {
		using type = std::false_type;
	};

	template<typename v1, typename v2, typename E = void>
	struct gt_helper {};

	template<typename v1, typename v2>
	struct gt_helper<v1, v2, std::enable_if_t<(v1::degree > v2::degree)>> {
		using type = std::true_type;
	};

	template<typename v1, typename v2>
	struct gt_helper<v1, v2, std::enable_if_t<(v1::degree == v2::degree)>> {
		using type = typename Ring::template gt_t<typename v1::aN, typename v2::aN>;
	};

	template<typename v1, typename v2>
	struct gt_helper<v1, v2, std::enable_if_t<(v1::degree < v2::degree)>> {
		using type = std::false_type;
	};

	template<typename coeff, size_t deg>
	struct monomial {
		using type = poly_mul_t<Ring, typename polynomial<Ring>::X, typename monomial<coeff, deg - 1>::type>;
	};

	template<typename coeff>
	struct monomial<coeff, 0> {
		using type = val<coeff>;
	};


public:
	using zero = typename polynomial<Ring>::template val<typename Ring::zero>;
	using one = typename polynomial<Ring>::template val<typename Ring::one>;
	using X = typename polynomial<Ring>::template val<typename Ring::one, typename Ring::zero>;

	template<typename v1, typename v2>
	using add_t = poly_add_t<Ring, v1, v2>;

	template<typename v1, typename v2>
	using sub_t = poly_sub_t<Ring, v1, v2>;

	template<typename v1, typename v2>
	using mul_t = poly_mul_t<Ring, v1, v2>;

	template<typename v1, typename v2>
	using eq_t = typename eq_helper<v1, v2>::type;

	template<typename v1, typename v2>
	using lt_t = typename lt_helper<v1, v2>::type;

	template<typename v1, typename v2>
	using gt_t = typename gt_helper<v1, v2>::type;

	template<typename v1, typename v2>
	using div_t = typename poly_div_helper<Ring, v1, v2, zero, v1>::type;

	template<typename coeff, size_t deg>
	using monomial_t = typename monomial<coeff, deg>::type;
};

template<typename Ring, typename P, typename E>
struct poly_simplify {};

// when high power is zero : strip
template<typename Ring, typename P>
struct poly_simplify<Ring, P, typename std::enable_if<
	std::is_same<
	typename Ring::zero,
	typename P::aN
	>::value && (P::degree > 0)
	>::type>
	{
		using type = typename poly_simplify<Ring, typename P::strip>::type;
	};

// otherwise : do nothing
template<typename Ring, typename P>
struct poly_simplify<Ring, P, typename std::enable_if<
	!std::is_same<
	typename Ring::zero,
	typename P::aN
	>::value && (P::degree > 0)
>::type>
{
	using type = P;
};

// do not simplify constants
template<typename Ring, typename P>
struct poly_simplify<Ring, P, std::enable_if_t<P::degree == 0>> {
	using type = P;
};

// division helper
template<typename Ring, typename A, typename B, typename Q, typename R, typename E>
struct poly_div_helper {};

template<typename Ring, typename A, typename B, typename Q, typename R>
struct poly_div_helper<Ring, A, B, Q, R, std::enable_if_t<
	(R::degree < B::degree) ||
	(R::degree == 0 && std::is_same<typename R::aN, typename Ring::zero>::value)>> {
	using type = Q;
};

template<typename Ring, typename A, typename B, typename Q, typename R>
struct poly_div_helper<Ring, A, B, Q, R, std::enable_if_t<
	(R::degree >= B::degree) &&
	!(R::degree == 0 && std::is_same<typename R::aN, typename Ring::zero>::value)>> {
private:
	using rN = typename R::aN;
	using bN = typename B::aN;
	using pT = typename polynomial<Ring>::template monomial_t<typename Ring::template div_t<rN, bN>, R::degree - B::degree>;

public:
	using type = typename poly_div_helper<
		Ring,
		A,
		B,
		typename polynomial<Ring>::template add_t<Q, pT>,
		typename polynomial<Ring>::template sub_t<R, typename polynomial<Ring>::template mul_t<pT, B>>>::type;
};

// addition at
template<typename Ring, typename P1, typename P2, size_t index>
struct poly_add_at {
	using type =
		typename Ring::template add_t<typename P1::template coeff_at_t<index>, typename P2::template coeff_at_t<index>>;
};

template<typename Ring, typename P1, typename P2, size_t index>
using poly_add_at_t = typename poly_add_at<Ring, P1, P2, index>::type;

template<typename Ring, typename P1, typename P2, std::size_t... I>
struct poly_add_low<Ring, P1, P2, std::index_sequence<I...>> {
	using type = typename polynomial<Ring>::template val<poly_add_at_t<Ring, P1, P2, I>...>;
};

// substraction at
template<typename Ring, typename P1, typename P2, size_t index>
struct poly_sub_at {
	using type =
		typename Ring::template sub_t<typename P1::template coeff_at_t<index>, typename P2::template coeff_at_t<index>>;
};

template<typename Ring, typename P1, typename P2, size_t index>
using poly_sub_at_t = typename poly_sub_at<Ring, P1, P2, index>::type;

template<typename Ring, typename P1, typename P2, std::size_t... I>
struct poly_sub_low<Ring, P1, P2, std::index_sequence<I...>> {
	using type = typename polynomial<Ring>::template val<poly_sub_at_t<Ring, P1, P2, I>...>;
};

// multiplication at
template<typename Ring, typename v1, typename v2, size_t k, size_t index, size_t stop>
struct poly_mul_at_loop_helper {
	using type = typename Ring::template add_t<
		typename Ring::template mul_t<
		typename v1::template coeff_at_t<index>,
		typename v2::template coeff_at_t<k - index>
		>,
		typename poly_mul_at_loop_helper<Ring, v1, v2, k, index + 1, stop>::type
	>;
};

template<typename Ring, typename v1, typename v2, size_t k, size_t stop>
struct poly_mul_at_loop_helper <Ring, v1, v2, k, stop, stop> {
	using type = typename Ring::template mul_t<typename v1::template coeff_at_t<stop>, typename v2::template coeff_at_t<0>>;
};

template <typename Ring, typename v1, typename v2, size_t k, typename E = void>
struct poly_mul_at {};

template<typename Ring, typename v1, typename v2, size_t k>
struct poly_mul_at<Ring, v1, v2, k, std::enable_if_t<(k < 0) || (k > v1::degree + v2::degree)>> {
	using type = Ring::zero;
};

template<typename Ring, typename v1, typename v2, size_t k>
struct poly_mul_at<Ring, v1, v2, k, std::enable_if_t<(k >= 0) && (k <= v1::degree + v2::degree)>> {
	using type = typename poly_mul_at_loop_helper<Ring, v1, v2, k, 0, k>::type;
};

template<typename Ring, typename P1, typename P2, size_t index>
using poly_mul_at_t = typename poly_mul_at<Ring, P1, P2, index>::type;

template<typename Ring, typename P1, typename P2, std::size_t... I>
struct poly_mul_low<Ring, P1, P2, std::index_sequence<I...>> {
	using type = typename polynomial<Ring>::template val<poly_mul_at_t<Ring, P1, P2, I>...>;
};

template<typename Ring>
struct _FractionField {
	static constexpr bool is_field = true;

	template<typename val1, typename val2>
	struct val {
		using x = val1;
		using y = val2;
	};

	using zero = val<typename Ring::zero, typename Ring::one>;
	using one = val<typename Ring::one, typename Ring::one>;

private:
	template<typename v, typename E = void>
	struct simplify {};

	template<typename v>
	struct simplify<v, std::enable_if_t<v::x::v != 0>> {
		using type = typename _FractionField<Ring>::template val<
			typename Ring::template div_t<
			typename v::x,
			typename Ring::template gcd_t<typename v::x, typename v::y>>,
			typename Ring::template div_t<
			typename v::y,
			typename Ring::template gcd_t<typename v::x, typename v::y>>>;
	};

	template<typename v>
	struct simplify<v, std::enable_if_t<v::x::v == 0>> {
		using type = typename _FractionField<Ring>::zero;
	};

	template<typename v>
	using simplify_t = simplify<v>::type;

	template<typename v1, typename v2>
	struct add {
	private:
		using a = typename Ring::template mul_t<typename v1::x, typename v2::y>;
		using b = typename Ring::template mul_t<typename v1::y, typename v2::x>;
		using dividend = typename Ring::template add_t<a, b>;
		using diviser = typename Ring::template mul_t<typename v1::y, typename v2::y>;
		using g = typename Ring::template gcd_t<dividend, diviser>;

	public:
		using type = typename _FractionField<Ring>::template simplify_t<val<dividend, diviser>>;
	};

	template<typename v1, typename v2>
	struct sub {
	private:
		using a = typename Ring::template mul_t<typename v1::x, typename v2::y>;
		using b = typename Ring::template mul_t<typename v1::y, typename v2::x>;
		using dividend = typename Ring::template sub_t<a, b>;
		using diviser = typename Ring::template mul_t<typename v1::y, typename v2::y>;
		using g = typename Ring::template gcd_t<dividend, diviser>;

	public:
		using type = typename _FractionField<Ring>::template simplify_t<val<dividend, diviser>>;
	};

	template<typename v1, typename v2>
	struct mul {
	private:
		using a = typename Ring::template mul_t<typename v1::x, typename v2::x>;
		using b = typename Ring::template mul_t<typename v1::y, typename v2::y>;

	public:
		using type = _FractionField<Ring>::template simplify_t<val<a, b>>;
	};

	template<typename v1, typename v2>
	struct div {
	private:
		using a = typename Ring::template mul_t<typename v1::x, typename v2::y>;
		using b = typename Ring::template mul_t<typename v1::y, typename v2::x>;

	public:
		using type = _FractionField<Ring>::template simplify_t<val<a, b>>;
	};

public:
	template<typename v1, typename v2>
	using add_t = typename add<v1, v2>::type;
	template<typename v1, typename v2>
	using sub_t = typename sub<v1, v2>::type;
	template<typename v1, typename v2>
	using mul_t = typename mul<v1, v2>::type;
	template<typename v1, typename v2>
	using div_t = typename div<v1, v2>::type;
};

template<typename Ring, typename E = void>
struct FractionFieldImpl {};

// fraction field of a field is the field itself
template<typename Field>
struct FractionFieldImpl<Field, std::enable_if_t<Field::is_field>> {
	using type = Field;
};

// fraction field of a ring is the actual fraction field
template<typename Ring>
struct FractionFieldImpl<Ring, std::enable_if_t<!Ring::is_field>> {
	using type = _FractionField<Ring>;
};

template<typename Ring>
using FractionField = typename FractionFieldImpl<Ring>::type;

using Q32 = FractionField<i32>;