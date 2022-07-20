namespace aerobus {
    namespace internal {
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
				using ring_type = Ring;
				using field_type = _FractionField<Ring>;

			 	static constexpr bool is_integer = std::is_same<val2, typename Ring::one>::value;

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

			template<typename v>
			using inject_t = val<v, typename Ring::one>;

			template<auto x>
			using inject_constant_t = val<typename Ring::template inject_constant_t<x>, typename Ring::one>;

			template<typename v>
			using inject_ring_t = val<typename Ring::template inject_ring_t<v>, typename Ring::one>;

			using ring_type = Ring;

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

				using posx = std::conditional_t<!Ring::template pos_t<newy>::value, typename Ring::template sub_t<typename Ring::zero, newx>, newx>;
				using posy = std::conditional_t<!Ring::template pos_t<newy>::value, typename Ring::template sub_t<typename Ring::zero, newy>, newy>;
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

			template<typename v>
			struct pos {
				using type = std::conditional_t<
					(Ring::template pos_t<typename v::x>::value && Ring::template pos_t<typename v::y>::value) ||
					(!Ring::template pos_t<typename v::x>::value && !Ring::template pos_t<typename v::y>::value),
					std::true_type,
					std::false_type>;

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

			template<typename v1, typename v2, typename E = void>
			struct div {};
			
			template<typename v1, typename v2>
			struct div<v1, v2, std::enable_if_t<!std::is_same<v2, typename _FractionField<Ring>::zero>::value>>  {
			private:
				using a = typename Ring::template mul_t<typename v1::x, typename v2::y>;
				using b = typename Ring::template mul_t<typename v1::y, typename v2::x>;

			public:
				using type = typename _FractionField<Ring>::template simplify_t<val<a, b>>;
			};

			template<typename v1, typename v2>
			struct div<v1, v2, std::enable_if_t<
				std::is_same<zero, v1>::value && std::is_same<v2, zero>::value>>  {
				using type = one;
			};

			template<typename v1, typename v2>
			struct eq {
				using type = std::conditional_t<
						std::is_same<typename simplify_t<v1>::x, typename simplify_t<v2>::x>::value &&
						std::is_same<typename simplify_t<v1>::y, typename simplify_t<v2>::y>::value,
					std::true_type, 
					std::false_type>;
			};

			template<typename TL, typename E = void>
			struct vadd {};

			template<typename TL>
			struct vadd<TL, std::enable_if_t<(TL::length > 1)>> {
				using head = typename TL::pop_front::type;
				using tail = typename TL::pop_front::tail;
				using type = typename add<head, typename vadd<tail>::type>::type;
			};

			template<typename TL>
			struct vadd<TL, std::enable_if_t<(TL::length == 1)>> {
				using type = typename TL::template at<0>;
			};

			template<typename... vals>
			struct vmul {};

			template<typename v1, typename... vals>
			struct vmul<v1, vals...> {
				using type = typename mul<v1, typename vmul<vals...>::type>::type;
			};

			template<typename v1>
			struct vmul<v1> {
				using type = v1;
			};
		

		public:
			template<typename v1, typename v2>
			using add_t = typename add<v1, v2>::type;
			template<typename... vs>
			using vadd_t = typename vadd<vs...>::type;
			template<typename... vs>
			using vmul_t = typename vmul<vs...>::type;
			template<typename v1, typename v2>
			using sub_t = typename sub<v1, v2>::type;
			template<typename v1, typename v2>
			using mul_t = typename mul<v1, v2>::type;
			template<typename v1, typename v2>
			using div_t = typename div<v1, v2>::type;
			template<typename v1, typename v2>
			using eq_t = typename eq<v1, v2>::type;
			template<typename v1>
			using pos_t = typename pos<v1>::type;
		};

		template<typename Ring, typename E = void>
		struct FractionFieldImpl {};

		// fraction field of a field is the field itself
		template<typename Field>
		struct FractionFieldImpl<Field, std::enable_if_t<Field::is_field>> {
			using type = Field;
			template<typename v>
			using inject_t = v;
		};

		// fraction field of a ring is the actual fraction field
		template<typename Ring>
		struct FractionFieldImpl<Ring, std::enable_if_t<!Ring::is_field>> {
			using type = _FractionField<Ring>;
		};
    }

	template<typename Ring>
	using FractionField = typename internal::FractionFieldImpl<Ring>::type;

	template<typename Ring, char variable_name = 'x'>
	struct RationalFraction: FractionField<polynomial<Ring, variable_name>> {
		static_assert(Ring::is_field, "invalid ring for rational fractions");

		private:
		template<typename v>
		struct remove_divisor {
			using type = std::conditional_t<v::y::degree == 0, 
				typename FractionField<polynomial<Ring, variable_name>>::template val<
					typename polynomial<Ring>::template div_t<
						typename v::x, 
						typename polynomial<Ring, variable_name>::template val<typename v::y::aN>
					>,
					typename polynomial<Ring, variable_name>::one
				>,
				v
			>;
		};

		public:
		template<typename v>
		using remove_divisor_t = typename remove_divisor<v>::type;
	};
}