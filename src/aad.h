#pragma once



namespace aerobus {
    using P64X = polynomial<Q64, 'x'>;
    using FPQ64X = RationalFraction<Q64, 'x'>;
    using P64T = polynomial<P64X, 't'>;

    namespace aad {
        template <typename derived_t, typename E = void> 
        struct base_traits;

        template<typename E>
        using simplify_t = typename base_traits<E>::type;

        template<typename E1, typename E2, typename E = void>
        struct eq_helper : std::false_type {};

        template<typename E1, typename E2>
        struct eq_helper<E1, E2, 
            std::enable_if_t<
                std::is_same<simplify_t<E1>, simplify_t<E2>>::value>>: std::true_type {};

        template<typename E1, typename E2>
        using eq_t = typename eq_helper<E1, E2>::type;

        template<typename Derived>
        struct Expression {
            using type = typename base_traits<Derived>::type;
            using value_at_t0 = typename base_traits<Derived>::value_at_t0;
            using derive_t = typename base_traits<Derived>::derive_t;
            using derive_x = typename base_traits<Derived>::derive_x;
            static constexpr uint64_t op_count = base_traits<Derived>::type::_op_count;

            static std::string to_string() {
                return type::_to_string();
            }
        };
    
        // val is a Q64 type
        template<typename val>
        struct ConstantExpression: Expression<ConstantExpression<val>> {
            template<typename>
            friend struct Expression;
            using v = val;
            private:
                static constexpr uint64_t _op_count = 1;
                static std::string _to_string() {
                    return "(" + val::to_string() + ")";
                }
        };

        template<typename E>
        struct MinExpression: Expression<MinExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;

            private:
                static constexpr uint64_t _op_count = 1;
                static std::string _to_string() {
                    return "-(" + E::to_string() + ")";
                }
        };

        // a*X^k
        template<typename a, size_t k>
        struct XExpression: Expression<XExpression<a, k>> {
            template<typename>
            friend struct Expression;

            private:
                static constexpr uint64_t _op_count = 1;
                static std::string _to_string() {
                    std::string result = "";
                    if (!a::is_zero_t::value && !std::is_same<a, typename Q64::one>::value) {
                        result += a::to_string();
                    }
                    result += "X";
                    if (k > 1) {
                        result += "^" + std::to_string(k);
                    }
                    return result;
                }
        };

        // a*t^k
        template<typename a, size_t k>
        struct TExpression: Expression<TExpression<a, k>> {
            template<typename>
            friend struct Expression;

            private:
                static constexpr uint64_t _op_count = 1;
                static std::string _to_string() {
                    std::string result = "";
                    if (!a::is_zero_t::value && !std::is_same<a, typename Q64::one>::value) {
                        result += a::to_string();
                    }
                    result += "t";
                    if (k > 1) {
                        result += "^" + std::to_string(k);
                    }
                    return result;
                }
        };

        template<char op, typename TL>
        struct to_string_helper {
            static std::string func() {
                std::string result = "";
                if (TL::length > 0)
                {   
                    result += TL::template at<0>::to_string();
                    if (TL::length > 1)
                    {
                        result += std::string(1, op) + to_string_helper<op, typename TL::pop_front::tail>::func();
                    }
                }

                return result;
            }
        };

        template<typename operand_list>
        struct AddExpression : Expression<AddExpression<operand_list>> {
            template<typename>
            friend struct Expression;
            
            private:
                static constexpr uint64_t _op_count = operand_list::length;
                static std::string _to_string() {
                    return to_string_helper<'+', operand_list>::func();
                }
        };

        template<typename E1, typename E2>
        struct SubExpression : Expression<SubExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            using left = E1;
            using right = E2;
            private:
                static constexpr uint64_t _op_count = 2;
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") - (" + E2::to_string() + "))";
                }
        };
        
        template<typename operand_list>
        struct MulExpression : Expression<MulExpression<operand_list>> {
            template<typename>
            friend struct Expression;

            static constexpr size_t op_count = operand_list::length;
            
            using operands = operand_list;

            template<size_t i>
            using at = typename operand_list::template at<i>;

            private:
                static constexpr uint64_t _op_count = operand_list::length;
                static std::string _to_string() {
                    return to_string_helper<'*', operand_list>::func();
                }
        };

        template<typename E1, typename E2>
        struct DivExpression : Expression<DivExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            using left = E1;
            using right = E2;
            private:
                static constexpr uint64_t _op_count = 2;
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") / (" + E2::to_string() + "))";
                }
        };

        // can't exist by itself but by composition
        template<typename E>
        struct ExpExpression : Expression<ExpExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;
            private:
                static constexpr uint64_t _op_count = 1;
                static std::string _to_string() {
                    return "exp(" + E::to_string() + ")";
                }
        }; 
        
        template<typename E>
        struct SinExpression : Expression<SinExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;
            private:
                static constexpr uint64_t _op_count = 1;
                static std::string _to_string() {
                    return "sin(" + E::to_string() + ")";
                }
        };
        
        template<typename E>
        struct CosExpression : Expression<CosExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;
            private:
                static constexpr uint64_t _op_count = 1;
                static std::string _to_string() {
                    return "cos(" + E::to_string() + ")";
                }
        };

        template<typename E1, typename E2>
        struct base_traits<DivExpression<E1, E2>> {
            using type = DivExpression<E1, E2>;
            using value_at_t0 = typename FPQ64X::template div_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = DivExpression<
                                SubExpression<
                                    MulExpression<type_list<typename base_traits<E1>::derive_t, E2>>,
                                    MulExpression<type_list<E1, typename base_traits<E2>::derive_t>>
                                >,
                                MulExpression<type_list<E2, E2>>
                            >;
            
            using derive_x = DivExpression<
                                SubExpression<
                                    MulExpression<type_list<typename base_traits<E1>::derive_x, E2>>,
                                    MulExpression<type_list<E1, typename base_traits<E2>::derive_x>>
                                >,
                                MulExpression<type_list<E2, E2>>
                            >;
        };

        template<typename val>
        struct base_traits<ConstantExpression<val>> {
            using type = ConstantExpression<val>;
            using value_at_t0 = typename FPQ64X::inject_ring_t<val>;
            using derive_x = ConstantExpression<val>;
            using derive_t = ConstantExpression<typename Q64::zero>;
        };

        template<typename a, size_t k>
        struct base_traits<XExpression<a, k>> {
            public:
                using type = XExpression<a, k>;
                using value_at_t0 =  FPQ64X::template val<polynomial<Q64>::monomial_t<a, k>, typename polynomial<Q64>::one>;
                using derive_t = ConstantExpression<typename Q64::zero>;
                using derive_x = std::conditional_t<(k == 0), ConstantExpression<typename Q64::zero>,
                                    std::conditional_t<k == 1, ConstantExpression<a>, 
                                        XExpression<Q64::mul_t<a, Q64::inject_constant_t<k>>, k-1>>>;
        };
        
        template<typename a, size_t k>
        struct base_traits<TExpression<a, k>> {
            using type = TExpression<a, k>;
            using value_at_t0 = std::conditional_t<
                k == 0,
                typename FPQ64X::inject_ring_t<a>,
                typename FPQ64X::zero
            >;
            using derive_t = std::conditional_t<(k == 0), 
                                ConstantExpression<typename Q64::zero>,
                                std::conditional_t<k == 1, 
                                        ConstantExpression<a>, 
                                        TExpression<Q64::mul_t<a, Q64::inject_constant_t<k>>, k-1>>>;
            using derive_x = ConstantExpression<typename Q64::zero>;
        };

        template<typename E>
        struct base_traits<ExpExpression<E>, 
            std::enable_if_t<!eq_t<ConstantExpression<typename Q64::zero>, simplify_t<E>>::value>
        > {
            using type = ExpExpression<E>;
            using value_at_t0 = typename FPQ64X::one; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<type_list<typename base_traits<E>::derive_t, ExpExpression<simplify_t<E>>>>;
            using derive_x = MulExpression<type_list<typename base_traits<E>::derive_x, ExpExpression<simplify_t<E>>>>;
        };
        
        template<typename E>
        struct base_traits<SinExpression<E>> {
            using type = SinExpression<E>;
            using value_at_t0 = typename FPQ64X::zero; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<type_list<typename base_traits<E>::derive_t, CosExpression<E>>>;
            using derive_x = MulExpression<type_list<typename base_traits<E>::derive_x, CosExpression<E>>>;
        };
        
        template<typename E>
        struct base_traits<CosExpression<E>> {
            using type = CosExpression<E>;
            using value_at_t0 = typename FPQ64X::one; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<type_list<typename base_traits<E>::derive_t, SubExpression<ConstantExpression<Q64::zero>, SinExpression<E>>>>;
            using derive_x = MulExpression<type_list<typename base_traits<E>::derive_x, SubExpression<ConstantExpression<Q64::zero>, SinExpression<E>>>>;
        };

        // export facilities
        using T = aad::TExpression<typename Q64::one, 1>;
        using X = aad::XExpression<typename  Q64::one, 1>;
        using XT = aad::MulExpression<type_list<X, T>>;
        using ONE = aad::ConstantExpression<typename Q64::one>;
        using ZERO = aad::ConstantExpression<typename Q64::zero>;
        using TWO = aad::ConstantExpression<Q64::inject_constant_t<2>>;

        // x - y
        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>> {
            using type = SubExpression<E1, E2>;
            using value_at_t0 = FPQ64X::template sub_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = SubExpression<typename base_traits<E1>::derive_t, typename base_traits<E2>::derive_t>;
            using derive_x = SubExpression<typename base_traits<E1>::derive_x, typename base_traits<E2>::derive_x>;
        };

        // -x
        template<typename E>
        struct base_traits<MinExpression<E>> {
            using type = MinExpression<E>;
            using value_at_t0 = FPQ64X::sub_t<typename FPQ64X::zero, typename base_traits<E>::value_at_t0>;
            using derive_t = MinExpression<typename base_traits<E>::derive_t>;
            using derive_x = MinExpression<typename base_traits<E>::derive_x>;
        };

        // // x * y
        // template<typename Op1, typename Op2, typename Op3, typename... Ops>
        // struct base_traits<MulExpression<Op1, Op2, Op3, Ops...>>
        // {
        //     using type = MulExpression<Op1, Op2, Op3, Ops...>;
        //     using value_at_t0 = FPQ64X::template vmul_t<
        //                             typename base_traits<Op1>::value_at_t0, 
        //                             typename base_traits<Op2>::value_at_t0, 
        //                             typename base_traits<Op3>::value_at_t0, 
        //                             typename base_traits<Ops>::value_at_t0...>;
        //     using derive_t = typename base_traits<MulExpression<Op1>, MulExpression<Op2, Op3, Ops...>>::derive_t;
        //     using derive_x = typename base_traits<MulExpression<Op1>, MulExpression<Op2, Op3, Ops...>>::derive_x;
        // };

        // template<typename E1, typename E2>
        // struct base_traits<MulExpression<E1, E2>> {
        //     using type = MulExpression<E1, E2>;
        //     using value_at_t0 = FPQ64X::template mul_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
        //     using derive_t = AddExpression<
        //                     MulExpression<typename base_traits<E1>::derive_t, E2>,
        //                     MulExpression<E1, typename base_traits<E2>::derive_t>>;
            
        //     using derive_x = AddExpression<
        //                     MulExpression<typename base_traits<E1>::derive_x, E2>,
        //                     MulExpression<E1, typename base_traits<E2>::derive_x>>;
        // };

        template<typename E, size_t k>
        struct taylor_coeff_helper {
            using type = FPQ64X::remove_divisor_t<FPQ64X::template div_t<typename taylor_coeff_helper<typename simplify_t<E>::derive_t, k-1>::type,
                                                FPQ64X::inject_constant_t<k>>>;
        };

        // binary add simplification
        template<typename E1, typename E2, typename E = void>
        struct simplify_binary_add;

        template<typename E1, typename E2>
        struct simplify_binary_add<E1, E2, 
                std::enable_if_t<
                    // two constants
                     (internal::is_instantiation_of<ConstantExpression, simplify_t<E1>>::value && 
                      internal::is_instantiation_of<ConstantExpression, simplify_t<E2>>::value)
                >
            > {
            using type = ConstantExpression<typename Q64::add_t<typename simplify_t<E1>::v, typename simplify_t<E2>::v>>;
        };

        template<typename E1, typename E2>
        struct simplify_binary_add<E1, E2, 
                std::enable_if_t<
                     // one zero, other not constant
                     !internal::is_instantiation_of<ConstantExpression, simplify_t<E2>>::value &&
                     std::is_same<ZERO, simplify_t<E1>>::value
                >
            > {
            using type = simplify_t<E2>;
        };

        template<typename E1, typename E2>
        struct simplify_binary_add<E1, E2, 
                std::enable_if_t<
                     // other zero, first not constant
                     !internal::is_instantiation_of<ConstantExpression, simplify_t<E1>>::value &&
                     std::is_same<ZERO, simplify_t<E2>>::value
                >
            > {
            using type = simplify_t<E1>;
        };

        template<typename E1, typename E2>
        struct simplify_binary_add<E1, E2, 
                std::enable_if_t<
                    // x + -x or -x + x
                    ((internal::is_instantiation_of<MinExpression, simplify_t<E1>>::value && std::is_same_v<typename simplify_t<E1>::v, simplify_t<E2>>) ||
                    (internal::is_instantiation_of<MinExpression, simplify_t<E2>>::value && std::is_same_v<typename simplify_t<E2>::v, simplify_t<E1>>))
                >
            > {
            using type = ZERO;
        };

        // template<>
        // struct simplify_binary_add<ONE, X> {
        //     using type = AddExpression<type_list<ONE, X>>;
        // };

        template<typename E1, typename E2>
        struct simplify_binary_add<E1, E2, 
                std::enable_if_t<
                    // general
                     !(
                        internal::is_instantiation_of_v<ConstantExpression, simplify_t<E1>> && 
                        internal::is_instantiation_of_v<ConstantExpression, simplify_t<E2>>
                       ) &&
                     !std::is_same_v<ZERO, simplify_t<E1>> &&
                     !std::is_same_v<ZERO, simplify_t<E2>> &&
                     !((internal::is_instantiation_of_v<MinExpression, simplify_t<E1>> && 
                        std::is_same_v<typename simplify_t<E1>::v, simplify_t<E2>>) ||
                      (internal::is_instantiation_of_v<MinExpression, simplify_t<E2>> && 
                        std::is_same_v<typename simplify_t<E2>::v, simplify_t<E1>>))
                >
            > {
            using type = AddExpression<type_list<E1, E2>>;
        };

        template<template<typename> typename AddE, typename TL, uint64_t op_count>
        struct simplify_add_h {
            template<uint64_t i, uint64_t j>
            struct inner {
                using si = typename TL::template at<i>;
                using sj = typename TL::template at<j>;
                using ss = typename simplify_binary_add<si, sj>::type;
                using type = std::conditional_t<
                    (ss::op_count < 2),
                    // then remove i and j and replace i by simplification
                    typename simplify_add_h<AddE, typename TL::template remove<j>::template remove<i>::template insert<i, ss>, op_count-1>::template inner<i, j>::type,
                    typename inner<i, j+1>::type
                >;
            };

            template<uint64_t i>
            struct inner<i, op_count> {
                using type = typename inner<i+1, i+2>::type;
            };

            template<uint64_t j>
            struct inner<op_count-1, j> {
                using type = AddE<TL>;
            };

            using type = typename inner<0, 1>::type;
        };

        template<template<typename> typename AddE, typename TL>
        struct simplify_add_h<AddE, TL, 2> {
            using type = typename simplify_binary_add<typename TL::template at<0>, typename TL::template at<1>>::type;
        };

        template<template<typename> typename AddE, typename TL>
        struct simplify_add_h<AddE, TL, 1> {
            using type = AddE<TL>;
        };
        
        template<template<typename> typename AddE, typename TL>
        struct simplify_add_h<AddE, TL, 0> {
            using type = AddE<TL>;
        };


        template<typename TL, uint64_t... Is>
        auto make_tl_value_at_t0(std::integer_sequence<uint64_t, Is...>) -> decltype(type_list<typename TL::template at<Is>::value_at_t0...> {});

        template<typename TL>
        using tl_value_at_t0 = decltype(make_tl_value_at_t0<TL>(std::make_integer_sequence<uint64_t, TL::length> {}));

        
        template<typename TL, uint64_t... Is>
        auto make_tl_derive_t(std::integer_sequence<uint64_t, Is...>) -> decltype(type_list<typename TL::template at<Is>::derive_t...> {});

        template<typename TL>
        using tl_derive_t = decltype(make_tl_derive_t<TL>(std::make_integer_sequence<uint64_t, TL::length> {}));
        
        template<typename TL, uint64_t... Is>
        auto make_tl_derive_x(std::integer_sequence<uint64_t, Is...>) -> decltype(type_list<typename TL::template at<Is>::derive_x...> {});

        template<typename TL>
        using tl_derive_x = decltype(make_tl_derive_x<TL>(std::make_integer_sequence<uint64_t, TL::length> {}));

        template <template <typename> typename AddE, typename TL>
        struct base_traits_add_h
        {
            template <typename> friend struct Expression;
        public:
            using type = typename simplify_add_h<AddE, TL, TL::length>::type;
            using value_at_t0 = FPQ64X::template vadd_t<tl_value_at_t0<TL>>;
            using derive_t = typename simplify_add_h<AddE, tl_derive_t<TL>, TL::length>::type;
            using derive_x = typename simplify_add_h<AddE, tl_derive_x<TL>, TL::length>::type;
        };

        template<typename TL>
        struct base_traits<AddExpression<TL>> : base_traits_add_h<AddExpression, TL> {};

        template<typename E>
        struct taylor_coeff_helper<E, 0> {
            using type = typename simplify_t<E>::value_at_t0;
        };

        template<typename E, size_t k>
        using taylor_coeff_t =  typename taylor_coeff_helper<E, k>::type::x;


        template<typename E, typename I>
        struct taylor_low {};

        template<typename E, size_t... k>
        struct taylor_low<E, std::index_sequence<k...>> {
            using type = P64T::template val<taylor_coeff_t<E, k>...>;
        };

        template<typename E, size_t k>
        using taylor_expansion_t = typename taylor_low<E, internal::make_index_sequence_reverse<k>>::type;
    }
}