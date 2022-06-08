#include <cstdint>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <algorithm>
#include <functional>
#include <string>

namespace aerobus {

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


	template<int32_t n, int32_t i, typename E = void>
	struct _is_prime {};

	template<int32_t n, int32_t i>
	struct _is_prime<n, i, std::enable_if_t<(n == 2 || n == 3)>> : std::true_type {};

	template<int32_t n, int32_t i>
	struct _is_prime<n, i, std::enable_if_t<(n != 2 && n != 3) && (n <= 1 || n % 2 == 0 || n % 3 == 0)>> : std::false_type {};

	struct _is_prime<n, i, std::enable_if_t<
		(n > 1 && n % 2 != 0 && n % 3 != 0) &&
		(i >= 5 && i * i <= n) &&
		(n% i == 0 || n % (i + 2) == 0)>> : std::false_type {};


	template<int32_t n, int32_t i>
	struct _is_prime<n, i, std::enable_if_t<
		(n > 1 && n % 2 != 0 && n % 3 != 0) &&
		(i >= 5 && i * i <= n) &&
		(n% i != 0 && n % (i + 2) != 0)>> {
		static constexpr bool value = _is_prime<n, i + 6>::value;
	};

	template<int32_t n, int32_t i>
	struct _is_prime<n, i, std::enable_if_t<
		(n > 1 && n % 2 != 0 && n % 3 != 0) &&
		(i >= 5 && i * i > n)>> : std::true_type {};


	template<int32_t n>
	struct is_prime {
		static constexpr bool value = _is_prime<n, 5>::value;
	};

	template <std::size_t ... Is>
	constexpr auto index_sequence_reverse(std::index_sequence<Is...> const&)
		-> decltype(std::index_sequence<sizeof...(Is) - 1U - Is...>{});

	template <std::size_t N>
	using make_index_sequence_reverse
		= decltype(index_sequence_reverse(std::make_index_sequence<N>{}));

	template<typename Ring, typename E = void>
	struct gcd;

	template<typename Ring>
	struct gcd<Ring, std::enable_if_t<Ring::is_integral_domain>> {
		template<typename A, typename B, typename E = void>
		struct gcd_helper {};

		// B = 0, A > 0
		template<typename A, typename B>
		struct gcd_helper<A, B, std::enable_if_t<
			((B::is_zero_t::value) &&
				(Ring::template gt_t<A, typename Ring::zero>::value))>>
		{
			using type = A;
		};

		// B = 0, A < 0
		template<typename A, typename B>
		struct gcd_helper<A, B, std::enable_if_t<
			((B::is_zero_t::value) &&
				!(Ring::template gt_t<A, typename Ring::zero>::value))>>
		{
			using type = typename Ring::template sub_t<typename Ring::zero, A>;
		};

		// B != 0
		template<typename A, typename B>
		struct gcd_helper<A, B, std::enable_if_t<
			(!B::is_zero_t::value)
			>> {
		private:
			// A / B
			using k = typename Ring::template div_t<A, B>;
			// A - (A/B)*B = A % B
			using m = typename Ring::template sub_t<A, typename Ring::template mul_t<k, B>>;
		public:
			using type = typename gcd_helper<B, m>::type;
		};

		template<typename A, typename B>
		using type = typename gcd_helper<A, B>::type;
	};

	struct i32 {
		using inner_type = int32_t;
		template<int32_t x>
		struct val {
			static constexpr int32_t v = x;

			template<typename valueType>
			static constexpr valueType get() { return static_cast<valueType>(x); }

			using is_zero_t = std::bool_constant<x == 0>;
			static std::string to_string() {
				return std::to_string(x);
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& v) {
				return static_cast<valueRing>(x);
			}
		};

		using zero = val<0>;
		using one = val<1>;
		static constexpr bool is_field = false;
		static constexpr bool is_integral_domain = true;

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
		struct remainder {
			using type = val<v1::v% v2::v>;
		};

		template<typename v1, typename v2>
		struct gt {
			using type = std::conditional_t<(v1::v > v2::v), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct lt {
			using type = std::conditional_t<(v1::v < v2::v), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct eq {
			using type = std::conditional_t<(v1::v == v2::v), std::true_type, std::false_type>;
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
		using mod_t = typename remainder<v1, v2>::type;

		template<typename v1, typename v2>
		using gt_t = typename gt<v1, v2>::type;

		template<typename v1, typename v2>
		using lt_t = typename lt<v1, v2>::type;

		template<typename v1, typename v2>
		using eq_t = typename eq<v1, v2>::type;

		template<typename v1, typename v2>
		using gcd_t = typename gcd<i32>::template type<v1, v2>;
	};

	struct i64 {
		using inner_type = int64_t;
		template<int64_t x>
		struct val {
			static constexpr int64_t v = x;

			template<typename valueType>
			static constexpr valueType get() { return static_cast<valueType>(x); }

			using is_zero_t = std::bool_constant<x == 0>;
			static std::string to_string() {
				return std::to_string(x);
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& v) {
				return static_cast<valueRing>(x);
			}
		};

		using zero = val<0>;
		using one = val<1>;
		static constexpr bool is_field = false;
		static constexpr bool is_integral_domain = true;

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
		struct remainder {
			using type = val<v1::v% v2::v>;
		};

		template<typename v1, typename v2>
		struct gt {
			using type = std::conditional_t<(v1::v > v2::v), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct lt {
			using type = std::conditional_t<(v1::v < v2::v), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct eq {
			using type = std::conditional_t<(v1::v == v2::v), std::true_type, std::false_type>;
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
		using mod_t = typename remainder<v1, v2>::type;

		template<typename v1, typename v2>
		using gt_t = typename gt<v1, v2>::type;

		template<typename v1, typename v2>
		using lt_t = typename lt<v1, v2>::type;

		template<typename v1, typename v2>
		using eq_t = typename eq<v1, v2>::type;

		template<typename v1, typename v2>
		using gcd_t = typename gcd<i64>::template type<v1, v2>;
	};

	template<int32_t p>
	struct zpz {
		using inner_type = int32_t;
		template<int32_t x>
		struct val {
			static constexpr int32_t v = x % p;

			template<typename valueType>
			static constexpr valueType get() { return static_cast<valueType>(x % p); }

			using is_zero_t = std::bool_constant<x% p == 0>;
			static std::string to_string() {
				return std::to_string(x % p);
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& v) {
				return static_cast<valueRing>(x % p);
			}
		};

		using zero = val<0>;
		using one = val<1>;
		static constexpr bool is_field = is_prime<p>::value;
		static constexpr bool is_integral_domain = is_field;

	private:
		template<typename v1, typename v2>
		struct add {
			using type = val<(v1::v + v2::v) % p>;
		};

		template<typename v1, typename v2>
		struct sub {
			using type = val<(v1::v - v2::v) % p>;
		};

		template<typename v1, typename v2>
		struct mul {
			using type = val<(v1::v* v2::v) % p>;
		};

		template<typename v1, typename v2>
		struct div {
			using type = val<(v1::v% p) / (v2::v % p)>;
		};

		template<typename v1, typename v2>
		struct remainder {
			using type = val<(v1::v% v2::v) % p>;
		};

		template<typename v1, typename v2>
		struct gt {
			using type = std::conditional_t<(v1::v% p > v2::v% p), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct lt {
			using type = std::conditional_t<(v1::v% p < v2::v% p), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct eq {
			using type = std::conditional_t<(v1::v% p == v2::v % p), std::true_type, std::false_type>;
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
		using mod_t = typename remainder<v1, v2>::type;

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

	template<typename Ring, typename coeff, typename... coeffs>
	struct poly_string_helper;

	template<typename coeffRing, typename valueRing, typename P>
	struct poly_eval_helper;

	template<typename Ring, typename P>
	struct make_unit;

	// coeffN x^N + ...
	template<typename Ring>
	struct polynomial {
		static constexpr bool is_field = false;
		static constexpr bool is_integral_domain = Ring::is_integral_domain;

		template<typename coeffN, typename... coeffs>
		struct val {
			static constexpr size_t degree = sizeof...(coeffs);
			using aN = coeffN;
			using strip = polynomial<Ring>::template val<coeffs...>;
			using is_zero_t = std::bool_constant<(degree == 0) && (aN::is_zero_t::value)>;

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

			static std::string to_string() {
				return poly_string_helper<Ring, coeffN, coeffs...>::func();
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& x) {
				return poly_eval_helper<Ring, valueRing, val>::template inner<0, degree + 1>::func(static_cast<valueRing>(0), x);
			}
		};

		// specialization for constants
		template<typename coeffN>
		struct val<coeffN> {
			static constexpr size_t degree = 0;
			using aN = coeffN;
			using strip = val<coeffN>;
			using is_zero_t = std::bool_constant<aN::is_zero_t::value>;

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

			static std::string to_string() {
				return poly_string_helper<Ring, coeffN>::apply();
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& x) {
				static_cast<valueRing>(aN::v);
			}
		};

		using zero = typename polynomial<Ring>::template val<typename Ring::zero>;
		using one = typename polynomial<Ring>::template val<typename Ring::one>;
		using X = typename polynomial<Ring>::template val<typename Ring::one, typename Ring::zero>;

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
			using type = std::false_type;
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
		using div_t = std::conditional_t<Ring::is_integral_domain, typename poly_div_helper<Ring, v1, v2, zero, v1>::q_type, void>;

		template<typename v1, typename v2>
		using mod_t = typename poly_div_helper<Ring, v1, v2, zero, v1>::mod_type;

		template<typename coeff, size_t deg>
		using monomial_t = typename monomial<coeff, deg>::type;

		template<typename v1, typename v2>
		using gcd_t = std::conditional_t<
			Ring::is_integral_domain,
			typename make_unit<Ring, typename gcd<polynomial<Ring>>::template type<v1, v2>>::type,
			void>;
	};

	template<typename coeffRing, typename valueRing, typename P>
	struct poly_eval_helper
	{
		template<size_t index, size_t stop>
		struct inner {
			static constexpr valueRing func(const valueRing& accum, const valueRing& x) {
				constexpr valueRing coeff = static_cast<valueRing>(P::template coeff_at_t<P::degree - index>::template get<valueRing>());
				return poly_eval_helper<coeffRing, valueRing, P>::template inner<index + 1, stop>::func(x * accum + coeff, x);
			}
		};

		template<size_t stop>
		struct inner<stop, stop> {
			static constexpr valueRing func(const valueRing& accum, const valueRing& x) {
				return accum;
			}
		};
	};


	template<typename Ring, typename P>
	struct make_unit {
		using type = typename polynomial<Ring>::template div_t<P, typename polynomial<Ring>::template val<typename P::aN>>;
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

	// to string helper -- use for debugging
	template<typename Ring, typename coeff, typename... coeffs>
	struct poly_string_helper {
		static std::string func() {
			if (Ring::template eq_t<coeff, typename Ring::one>::value) {
				if (sizeof...(coeffs) == 1) {
					return "x + " + poly_string_helper<Ring, coeffs...>::func();
				}
				else {
					return "x^" + std::to_string(sizeof...(coeffs)) + " + " + poly_string_helper<Ring, coeffs...>::func();
				}
			}
			else {
				if (sizeof...(coeffs) == 1) {
					return coeff::to_string() + " x" + " + " + poly_string_helper<Ring, coeffs...>::func();
				}
				else {
					return coeff::to_string() + " x^" + std::to_string(sizeof...(coeffs)) + " + " + poly_string_helper<Ring, coeffs...>::func();
				}
			}
		}
	};

	template<typename Ring, typename coeff>
	struct poly_string_helper<Ring, coeff> {
		static std::string func() {
			return coeff::to_string();
		}
	};

	// division helper
	template<typename Ring, typename A, typename B, typename Q, typename R, typename E>
	struct poly_div_helper {};

	template<typename Ring, typename A, typename B, typename Q, typename R>
	struct poly_div_helper<Ring, A, B, Q, R, std::enable_if_t<
		(R::degree < B::degree) ||
		(R::degree == 0 && std::is_same<typename R::aN, typename Ring::zero>::value)>> {
		using q_type = Q;
		using mod_type = R;
		using gcd_type = B;
	};

	template<typename Ring, typename A, typename B, typename Q, typename R>
	struct poly_div_helper<Ring, A, B, Q, R, std::enable_if_t<
		(R::degree >= B::degree) &&
		!(R::degree == 0 && std::is_same<typename R::aN, typename Ring::zero>::value)>> {
	private:
		using rN = typename R::aN;
		using bN = typename B::aN;
		using pT = typename polynomial<Ring>::template monomial_t<typename Ring::template div_t<rN, bN>, R::degree - B::degree>;
		using rr = typename polynomial<Ring>::template sub_t<R, typename polynomial<Ring>::template mul_t<pT, B>>;
		using qq = typename polynomial<Ring>::template add_t<Q, pT>;

	public:
		using q_type = typename poly_div_helper<Ring, A, B, qq, rr>::q_type;
		using mod_type = typename poly_div_helper<Ring, A, B, qq, rr>::mod_type;
		using gcd_type = rr;
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
		using type = typename Ring::zero;
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

	template<typename Ring, typename E = void>
	struct _FractionField {};

	template<typename Ring>
	struct _FractionField<Ring, std::enable_if_t<Ring::is_integral_domain>>
	{
		static constexpr bool is_field = true;
		static constexpr bool is_integral_domain = true;

		template<typename val1, typename val2, typename E = void>
		struct to_string_helper {};

		template<typename val1, typename val2>
		struct to_string_helper <val1, val2,
			std::enable_if_t<
			Ring::template eq_t<
			val2, typename Ring::one
			>::value
			>
		> {
			static std::string func() {
				return val1::to_string();
			}
		};

		template<typename val1, typename val2>
		struct to_string_helper<val1, val2,
			std::enable_if_t<
			!Ring::template eq_t<
			val2,
			typename Ring::one
			>::value
			>
		> {
			static std::string func() {
				return "(" + val1::to_string() + ") / (" + val2::to_string() + ")";
			}
		};

		template<typename val1, typename val2>
		struct val {
			using x = val1;
			using y = val2;
			using is_zero_t = typename val1::is_zero_t;

			template<typename valueType>
			static constexpr valueType get() { return static_cast<valueType>(x::v) / static_cast<valueType>(y::v); }

			static std::string to_string() {
				return to_string_helper<val1, val2>::func();
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& v) {
				return x::eval(v) / y::eval(v);
			}
		};

		using zero = val<typename Ring::zero, typename Ring::one>;
		using one = val<typename Ring::one, typename Ring::one>;

	private:
		template<typename v, typename E = void>
		struct simplify {};

		// x = 0
		template<typename v>
		struct simplify<v, std::enable_if_t<v::x::is_zero_t::value>> {
			using type = typename _FractionField<Ring>::zero;
		};

		// x != 0
		template<typename v>
		struct simplify<v, std::enable_if_t<!v::x::is_zero_t::value>> {

		private:
			using _gcd = typename Ring::template gcd_t<typename v::x, typename v::y>;
			using newx = typename Ring::template div_t<typename v::x, _gcd>;
			using newy = typename Ring::template div_t<typename v::y, _gcd>;

			using posx = std::conditional_t<Ring::template lt_t<newy, typename Ring::zero>::value, typename Ring::template sub_t<typename Ring::zero, newx>, newx>;
			using posy = std::conditional_t<Ring::template lt_t<newy, typename Ring::zero>::value, typename Ring::template sub_t<typename Ring::zero, newy>, newy>;
		public:
			using type = typename _FractionField<Ring>::template val<posx, posy>;
		};

	public:

		template<typename v>
		using simplify_t = typename simplify<v>::type;

	private:

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
			using type = typename _FractionField<Ring>::template simplify_t<val<a, b>>;
		};

		template<typename v1, typename v2>
		struct div {
		private:
			using a = typename Ring::template mul_t<typename v1::x, typename v2::y>;
			using b = typename Ring::template mul_t<typename v1::y, typename v2::x>;

		public:
			using type = typename _FractionField<Ring>::template simplify_t<val<a, b>>;
		};

		template<typename v1, typename v2>
		struct eq {
			using type = std::conditional_t<(simplify_t<v1>::x::v == simplify_t<v2>::x::v) && (simplify_t<v1>::y::v == simplify_t<v2>::y::v), std::true_type, std::false_type>;
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
		using eq_t = typename eq<v1, v2>::type;
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
	using FPQ32 = FractionField<polynomial<Q32>>;


	template<typename T, auto x>
	struct factorial {
	private:
		template<typename, auto>
		friend struct factorial;
	public:
		using type = typename T::template mul_t<typename T::template val<x>, typename factorial<T, x - 1>::type>;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	template<typename T>
	struct factorial<T, 0> {
	private:
		template<typename, auto>
		friend struct factorial;
	public:
		using type = typename T::one;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	template<typename T, auto k, auto n, typename E = void>
	struct combination_helper {};

	template<typename T, auto k, auto n>
	struct combination_helper<T, k, n, std::enable_if_t<(n >= 0 && k <= (n / 2) && k > 0)>> {
		using type = typename FractionField<T>::template mul_t<
			typename combination_helper<T, k - 1, n - 1>::type,
			typename FractionField<T>::template val<typename T::template val<n>, typename T::template val<k>>>;
	};

	template<typename T, auto k, auto n>
	struct combination_helper<T, k, n, std::enable_if_t<(n >= 0 && k > (n / 2) && k > 0)>> {
		using type = typename combination_helper<T, n - k, n>::type;
	};

	template<typename T, auto n>
	struct combination_helper<T, 0, n> {
		using type = typename FractionField<T>::one;
	};

	template<typename T, auto k, auto n>
	struct combination {
		using type = typename combination_helper<T, k, n>::type::x;
		static constexpr typename T::inner_type value = combination_helper<T, k, n>::type::template get<typename T::inner_type>();
	};


	template<typename T, auto m>
	struct bernouilli;

	template<typename T, typename accum, auto k, auto m>
	struct bernouilli_helper {
		using type = typename bernouilli_helper<
			T,
			typename FractionField<T>::template add_t<
			accum,
			typename FractionField<T>::template mul_t<
			typename FractionField<T>::template val<
			typename combination<T, k, m + 1>::type,
			typename T::one>,
			typename bernouilli<T, k>::type
			>
			>,
			k + 1,
			m
		>::type;
	};

	template<typename T, typename accum, auto m>
	struct bernouilli_helper<T, accum, m, m>
	{
		using type = accum;
	};

	template<typename T, auto m>
	struct bernouilli {
		using type = typename FractionField<T>::template mul_t<
			typename bernouilli_helper<T, typename FractionField<T>::zero, 0, m>::type,
			typename FractionField<T>::template val<
			typename T::template val<static_cast<typename T::inner_type>(-1)>,
			typename T::template val<static_cast<typename T::inner_type>(m + 1)>
			>
		>;

		template<typename floatType>
		static constexpr floatType value = type::template get<floatType>();
	};

	template<typename T>
	struct bernouilli<T, 0> {
		using type = typename FractionField<T>::one;

		template<typename floatType>
		static constexpr floatType value = type::template get<floatType>();
	};

	template<typename T, int k, typename E = void>
	struct alternate {};

	template<typename T, int k>
	struct alternate<T, k, std::enable_if_t<k % 2 == 0>> {
		using type = typename T::one;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	template<typename T, int k>
	struct alternate<T, k, std::enable_if_t<k % 2 != 0>> {
		using type = typename T::template sub_t<T::zero, T::one>;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	template<typename T, auto p, auto n>
	struct pow { 
		using type = typename T::template mul_t<typename T::template val<p>, typename pow<T, p, n-1>::type>;
	};

	template<typename T, auto p>
	struct pow<T, p, 0> { using type = typename T::one; };


	template<typename, template<typename, size_t> typename, class>
	struct make_taylor_impl;

	template<typename T, template<typename, size_t> typename coeff_at, size_t... Is>
	struct make_taylor_impl<T, coeff_at, std::integer_sequence<size_t, Is...>> {
		using type = typename polynomial<FractionField<T>>::template val<typename coeff_at<T, Is>::type...>;
	};

	// generic taylor serie, depending on coefficients
	template<typename T, template<typename, size_t index> typename coeff_at, size_t deg>
	using taylor = typename make_taylor_impl<T, coeff_at, make_index_sequence_reverse<deg + 1>>::type;

	template<typename T, size_t i>
	struct exp_coeff {
		using type = typename FractionField<T>::template val<typename T::one, typename factorial<T, i>::type>;
	};

	template<typename T, size_t i, typename E = void>
	struct sin_coeff_helper {};

	template<typename T, size_t i>
	struct sin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct sin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = typename FractionField<T>::template val<typename alternate<T, i / 2>::type, typename factorial<T, i>::type>;
	};

	template<typename T, size_t i>
	struct sin_coeff {
		using type = typename sin_coeff_helper<T, i>::type;
	};

	template<typename T, size_t i, typename E = void>
	struct sh_coeff_helper {};

	template<typename T, size_t i>
	struct sh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct sh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = typename FractionField<T>::template val<typename T::one, typename factorial<T, i>::type>;
	};

	template<typename T, size_t i>
	struct sh_coeff {
		using type = typename sh_coeff_helper<T, i>::type;
	};

	template<typename T, size_t i, typename E = void>
	struct cos_coeff_helper {};

	template<typename T, size_t i>
	struct cos_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct cos_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = typename FractionField<T>::template val<typename alternate<T, i / 2>::type, typename factorial<T, i>::type>;
	};

	template<typename T, size_t i>
	struct cos_coeff {
		using type = typename cos_coeff_helper<T, i>::type;
	};

	template<typename T, size_t i, typename E = void>
	struct cosh_coeff_helper {};

	template<typename T, size_t i>
	struct cosh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct cosh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = typename FractionField<T>::template val<typename T::one, typename factorial<T, i>::type>;
	};

	template<typename T, size_t i>
	struct cosh_coeff {
		using type = typename cosh_coeff_helper<T, i>::type;
	};

	template<typename T, size_t i>
	struct geom_coeff { using type = typename FractionField<T>::one; };


	template<typename T, size_t i, typename E = void>
	struct atan_coeff_helper;

	template<typename T, size_t i>
	struct atan_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = typename FractionField<T>::template val<alternate<T, i / 2>::type, typename T::template val<i>>;
	};

	template<typename T, size_t i>
	struct atan_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct atan_coeff { using type = typename atan_coeff_helper<T, i>::type; };

	template<typename T, size_t i, typename E = void>
	struct asin_coeff_helper;

	template<typename T, size_t i>
	struct asin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type>
	{
		using type = typename FractionField<T>::template val<
			typename factorial<T, i - 1>::type,
			typename T::template mul_t<
				typename T::template val<i>,
					typename T::template mul_t<
						typename pow<T, 4, i / 2>::type,
						pow<T, factorial<T, i / 2>::value, 2>
			>>>;
	};

	template<typename T, size_t i>
	struct asin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type>
	{
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct asin_coeff { 
		using type = typename asin_coeff_helper<T, i>::type;
	};

	template<typename T, size_t i>
	struct lnp1_coeff { 
		using type = typename FractionField<T>::template val<
			typename alternate<T, i + 1>::type,
			typename T::template val<i>>;
	};

	template<typename T>
	struct lnp1_coeff<T, 0> { using type = typename FractionField<T>::zero; };

	template<typename T, size_t i, typename E = void>
	struct asinh_coeff_helper;

	template<typename T, size_t i>
	struct asinh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type>
	{
		using type = typename FractionField<T>::template val<
			typename T::template mul_t<
				typename alternate<T, i/2>::type,
				typename factorial<T, i-1>::type
			>,
			typename T::template mul_t<
				typename T::template mul_t<
					typename T::template val<i>,
					typename pow<T, (factorial<T, i/2>::value), 2>::type
				>,
				typename pow<T, 4, i/2>::type
			>
		>;
	};

	template<typename T, size_t i>
	struct asinh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type>
	{
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct asinh_coeff { 
		using type = typename asinh_coeff_helper<T, i>::type; 
	};

	template<typename T, size_t i, typename E = void>
	struct atanh_coeff_helper;

	template<typename T, size_t i>
	struct atanh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type>
	{
		// 1/i
		using type = typename FractionField<T>:: template val<
			typename T::one, 
			typename T::template val<static_cast<typename T::inner_type>(i)>>;
	};

	template<typename T, size_t i>
	struct atanh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type>
	{
		using type = typename FractionField<T>::zero;
	};

	template<typename T, size_t i>
	struct atanh_coeff {
		using type = typename asinh_coeff_helper<T, i>::type; 
	};

	// e^x
	template<typename T, size_t deg>
	using exp = taylor<T, exp_coeff, deg>;

	// e^x - 1
	template<typename T, size_t deg>
	using expm1 = typename polynomial<FractionField<T>>::template sub_t<
		exp<T, deg>,
		typename polynomial<FractionField<T>>::one>;

	/// ln(1+x)
	template<typename T, size_t deg>
	using lnp1 = taylor<T, lnp1_coeff, deg>;

	/// atan(x);
	template<typename T, size_t deg>
	using atan = taylor<T, atan_coeff, deg>;

	/// sin(x)
	template<typename T, size_t deg>
	using sin = taylor<T, sin_coeff, deg>;

	/// sh(x)
	template<typename T, size_t deg>
	using sinh = taylor<T, sh_coeff, deg>;

	/// ch(x)
	template<typename T, size_t deg>
	using cosh = taylor<T, cosh_coeff, deg>;

	/// cos(x)
	template<typename T, size_t deg>
	using cos = taylor<T, cos_coeff, deg>;

	/// 1 / (1-x)
	template<typename T, size_t deg>
	using geometric_sum = taylor<T, geom_coeff, deg>;

	/// asin(x)
	template<typename T, size_t deg>
	using asin = taylor<T, asin_coeff, deg>;

	/// asinh(x)
	template<typename T, size_t deg>
	using asinh = taylor<T, asinh_coeff, deg>;

	/// atanh(x)
	template<typename T, size_t deg>
	using atanh = taylor<T, atanh_coeff, deg>;
}