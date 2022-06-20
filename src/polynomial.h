namespace aerobus {
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
				using type = internal::type_at_t<sizeof...(coeffs) - index, coeffN, coeffs...>;
			};

			template<size_t index>
			struct coeff_at<index, std::enable_if_t<(index < 0 || index > sizeof...(coeffs))>> {
				using type = typename Ring::zero;
			};

			template<size_t index>
			using coeff_at_t = typename coeff_at<index>::type;

			static std::string to_string() {
				return string_helper<coeffN, coeffs...>::func();
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& x) {
				return eval_helper<valueRing, val>::template inner<0, degree + 1>::func(static_cast<valueRing>(0), x);
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
				return string_helper<coeffN>::func();
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& x) {
				return static_cast<valueRing>(aN::template get<valueRing>());
			}
		};

		using zero = typename polynomial<Ring>::template val<typename Ring::zero>;
		using one = typename polynomial<Ring>::template val<typename Ring::one>;
		using X = typename polynomial<Ring>::template val<typename Ring::one, typename Ring::zero>;

	private:
		template<typename P, typename E = void>
		struct simplify;

		template <typename P1, typename P2, typename I>
		struct add_low;

		template<typename P1, typename P2>
		struct add {
			using type = typename simplify<typename add_low<
			P1,
			P2,
			internal::make_index_sequence_reverse<
			std::max(P1::degree, P2::degree) + 1
			>>::type>::type;
		};

		template <typename P1, typename P2, typename I>
		struct sub_low;

		template <typename P1, typename P2, typename I>
		struct mul_low;

		template<typename v1, typename v2>
		struct mul {
				using type = typename mul_low<
					v1,
					v2,
					internal::make_index_sequence_reverse<
					v1::degree + v2::degree + 1
					>>::type;
		};

		template<typename coeff, size_t deg>
		struct monomial;

		template<typename v, typename E = void>
		struct derive_helper {};

		template<typename v>
		struct derive_helper<v, std::enable_if_t<v::degree == 0>> {
			using type = zero;
		};

		template<typename v>
		struct derive_helper<v, std::enable_if_t<v::degree != 0>> {
			using type = typename add<typename derive_helper<typename v::strip>::type,
				typename monomial<
					typename Ring::template mul_t<
						typename v::aN,
						typename Ring::template inject_constant_t<(v::degree)>
					>,
					v::degree - 1
				>::type
			>::type;
		};

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

		// when high power is zero : strip
		template<typename P>
		struct simplify<P, typename std::enable_if<
			std::is_same<
			typename Ring::zero,
			typename P::aN
			>::value && (P::degree > 0)
		>::type>
		{
			using type = typename simplify<typename P::strip>::type;
		};

		// otherwise : do nothing
		template<typename P>
		struct simplify<P, typename std::enable_if<
			!std::is_same<
			typename Ring::zero,
			typename P::aN
			>::value && (P::degree > 0)
		>::type>
		{
			using type = P;
		};

		// do not simplify constants
		template<typename P>
		struct simplify<P, std::enable_if_t<P::degree == 0>> {
			using type = P;
		};

		// addition at
		template<typename P1, typename P2, size_t index>
		struct add_at {
			using type =
				typename Ring::template add_t<typename P1::template coeff_at_t<index>, typename P2::template coeff_at_t<index>>;
		};

		template<typename P1, typename P2, size_t index>
		using add_at_t = typename add_at<P1, P2, index>::type;

		template<typename P1, typename P2, std::size_t... I>
		struct add_low<P1, P2, std::index_sequence<I...>> {
			using type = val<add_at_t<P1, P2, I>...>;
		};

		// substraction at
		template<typename P1, typename P2, size_t index>
		struct sub_at {
			using type =
				typename Ring::template sub_t<typename P1::template coeff_at_t<index>, typename P2::template coeff_at_t<index>>;
		};

		template<typename P1, typename P2, size_t index>
		using sub_at_t = typename sub_at<P1, P2, index>::type;

		template<typename P1, typename P2, std::size_t... I>
		struct sub_low<P1, P2, std::index_sequence<I...>> {
			using type = val<sub_at_t<P1, P2, I>...>;
		};

		template<typename P1, typename P2>
		struct sub {
			using type = typename simplify<typename sub_low<
			P1,
			P2,
			internal::make_index_sequence_reverse<
			std::max(P1::degree, P2::degree) + 1
			>>::type>::type;
		};

		// multiplication at
		template<typename v1, typename v2, size_t k, size_t index, size_t stop>
		struct mul_at_loop_helper {
			using type = typename Ring::template add_t<
				typename Ring::template mul_t<
				typename v1::template coeff_at_t<index>,
				typename v2::template coeff_at_t<k - index>
				>,
				typename mul_at_loop_helper<v1, v2, k, index + 1, stop>::type
			>;
		};

		template<typename v1, typename v2, size_t k, size_t stop>
		struct mul_at_loop_helper<v1, v2, k, stop, stop> {
			using type = typename Ring::template mul_t<typename v1::template coeff_at_t<stop>, typename v2::template coeff_at_t<0>>;
		};

		template <typename v1, typename v2, size_t k, typename E = void>
		struct mul_at {};

		template<typename v1, typename v2, size_t k>
		struct mul_at<v1, v2, k, std::enable_if_t<(k < 0) || (k > v1::degree + v2::degree)>> {
			using type = typename Ring::zero;
		};

		template<typename v1, typename v2, size_t k>
		struct mul_at<v1, v2, k, std::enable_if_t<(k >= 0) && (k <= v1::degree + v2::degree)>> {
			using type = typename mul_at_loop_helper<v1, v2, k, 0, k>::type;
		};

		template<typename P1, typename P2, size_t index>
		using mul_at_t = typename mul_at<P1, P2, index>::type;

		template<typename P1, typename P2, std::size_t... I>
		struct mul_low<P1, P2, std::index_sequence<I...>> {
			using type = val<mul_at_t<P1, P2, I>...>;
		};

		// division helper
		template< typename A, typename B, typename Q, typename R, typename E = void>
		struct div_helper {};

		template<typename A, typename B, typename Q, typename R>
		struct div_helper<A, B, Q, R, std::enable_if_t<
			(R::degree < B::degree) ||
			(R::degree == 0 && std::is_same<typename R::aN, typename Ring::zero>::value)>> {
			using q_type = Q;
			using mod_type = R;
			using gcd_type = B;
		};

		template<typename A, typename B, typename Q, typename R>
		struct div_helper<A, B, Q, R, std::enable_if_t<
			(R::degree >= B::degree) &&
			!(R::degree == 0 && std::is_same<typename R::aN, typename Ring::zero>::value)>> {
		private:
			using rN = typename R::aN;
			using bN = typename B::aN;
			using pT = typename monomial<typename Ring::template div_t<rN, bN>, R::degree - B::degree>::type;
			using rr = typename sub<R, typename mul<pT, B>::type>::type;
			using qq = typename add<Q, pT>::type;

		public:
			using q_type = typename div_helper<A, B, qq, rr>::q_type;
			using mod_type = typename div_helper<A, B, qq, rr>::mod_type;
			using gcd_type = rr;
		};

		template<typename A, typename B>
		struct div {
			static_assert(Ring::is_integral_domain, "cannot divide in that type of Ring");
			using q_type = typename div_helper<A, B, zero, A>::q_type;
			using m_type = typename div_helper<A, B, zero, A>::mod_type;
		};

		
		template<typename P>
		struct make_unit {
			using type = typename div<P, val<typename P::aN>>::q_type;
		};

		template<typename coeff, size_t deg>
		struct monomial {
			using type = typename mul<X, typename monomial<coeff, deg - 1>::type>::type;
		};

		template<typename coeff>
		struct monomial<coeff, 0> {
			using type = val<coeff>;
		};

		template<typename valueRing, typename P>
		struct eval_helper
		{
			template<size_t index, size_t stop>
			struct inner {
				static constexpr valueRing func(const valueRing& accum, const valueRing& x) {
					constexpr valueRing coeff = static_cast<valueRing>(P::template coeff_at_t<P::degree - index>::template get<valueRing>());
					return eval_helper<valueRing, P>::template inner<index + 1, stop>::func(x * accum + coeff, x);
				}
			};

			template<size_t stop>
			struct inner<stop, stop> {
				static constexpr valueRing func(const valueRing& accum, const valueRing& x) {
					return accum;
				}
			};
		};

		template<typename coeff, typename... coeffs>
		struct string_helper {
			static std::string func() {
				if (Ring::template eq_t<coeff, typename Ring::zero>::value) {
					if (sizeof...(coeffs) == 1) {
						return string_helper<coeffs...>::func();
					}
					else {
						return string_helper<coeffs...>::func();
					}
				}
				else if (Ring::template eq_t<coeff, typename Ring::one>::value) {
					if (sizeof...(coeffs) == 1) {
						return "x + " + string_helper<coeffs...>::func();
					}
					else {
						return "x^" + std::to_string(sizeof...(coeffs)) + " + " + string_helper<coeffs...>::func();
					}
				}
				else {
					if (sizeof...(coeffs) == 1) {
						return coeff::to_string() + " x" + " + " + string_helper<coeffs...>::func();
					}
					else {
						return coeff::to_string() + " x^" + std::to_string(sizeof...(coeffs)) + " + " + string_helper<coeffs...>::func();
					}
				}
			}
		};

		template<typename coeff>
		struct string_helper<coeff> {
			static std::string func() {
				return coeff::to_string();
			}
		};

	public:
		template<typename P>
		using simplify_t = typename simplify<P>::type;

		template<typename v1, typename v2>
		using add_t = typename add<v1, v2>::type;

		template<typename v1, typename v2>
		using sub_t = typename sub<v1, v2>::type;

		template<typename v1, typename v2>
		using mul_t = typename mul<v1, v2>::type;

		template<typename v1, typename v2>
		using eq_t = typename eq_helper<v1, v2>::type;

		template<typename v1, typename v2>
		using lt_t = typename lt_helper<v1, v2>::type;

		template<typename v1, typename v2>
		using gt_t = typename gt_helper<v1, v2>::type;

		template<typename v1, typename v2>
		using div_t = typename div<v1, v2>::q_type;

		template<typename v1, typename v2>
		using mod_t = typename div_helper<v1, v2, zero, v1>::mod_type;

		template<typename coeff, size_t deg>
		using monomial_t = typename monomial<coeff, deg>::type;

		template<typename v>
		using derive_t = typename derive_helper<v>::type;

		template<typename v>
		using pos_t = typename Ring::template pos_t<typename v::aN>;

		template<typename v1, typename v2>
		using gcd_t = std::conditional_t<
			Ring::is_integral_domain,
			typename make_unit<typename internal::gcd<polynomial<Ring>>::template type<v1, v2>>::type,
			void>;

		template<auto x>
		using inject_constant_t = val<typename Ring::template inject_constant_t<x>>;

		template<typename v>
		using inject_ring_t = val<v>;
	};
}