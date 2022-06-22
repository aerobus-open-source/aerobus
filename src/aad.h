#pragma once

namespace aerobus {
    using P64X = polynomial<Q64, 'x'>;
    using FPQ64X = RationalFraction<Q64, 'x'>;
    using P64T = polynomial<P64X, 't'>;

    namespace aad {
        template <typename derived_t, typename E = void> 
        struct base_traits;

        template<typename Derived>
        struct Expression {
            using type = typename base_traits<Derived>::type;
            using value_at_t0 = typename base_traits<Derived>::value_at_t0;
            using derive_t = typename base_traits<Derived>::derive_t;
            using derive_x = typename base_traits<Derived>::derive_x;

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
                static std::string _to_string() {
                    return "(" + val::to_string() + ")";
                }
        };

        // a*X^k
        template<typename a, size_t k>
        struct XExpression: Expression<XExpression<a, k>> {
            template<typename>
            friend struct Expression;

            using type = XExpression<a, k>;
            private:
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

        template<typename E1, typename E2>
        struct AddExpression : Expression<AddExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") + (" + E2::to_string() + "))";
                }
        };

        template<typename E1, typename E2>
        struct SubExpression : Expression<SubExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") - (" + E2::to_string() + "))";
                }
        };
        
        template<typename E1, typename E2>
        struct MulExpression : Expression<MulExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") * (" + E2::to_string() + "))";
                }
        };

        template<typename E1, typename E2>
        struct DivExpression : Expression<DivExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") / (" + E2::to_string() + "))";
                }
        };

        // can't exist by itself but by composition
        template<typename E>
        struct ExpExpression : Expression<ExpExpression<E>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "exp(" + E::to_string() + ")";
                }
        }; 
        
        template<typename E>
        struct SinExpression : Expression<SinExpression<E>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "sin(" + E::to_string() + ")";
                }
        };
        
        template<typename E>
        struct CosExpression : Expression<CosExpression<E>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "cos(" + E::to_string() + ")";
                }
        };

        // base_traits specialization and simplification rules
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

        // x + y
        // y != 0
        // x not constant
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>, 
                    std::enable_if_t<
                        !std::is_same<typename base_traits<E2>::type, ConstantExpression<Q64::zero>>::value &&
                        !internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value
        >> {
            using type = AddExpression<E1, E2>;
            using value_at_t0 = FPQ64X::template add_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = AddExpression<typename base_traits<E1>::derive_t, typename base_traits<E2>::derive_t>;
            using derive_x = AddExpression<typename base_traits<E1>::derive_x, typename base_traits<E2>::derive_x>;
        };

        // x+0= x
        // x not constant
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                    std::enable_if_t<
                        std::is_same<typename base_traits<E2>::type, ConstantExpression<Q64::zero>>::value && 
                        !internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value 
                >> {
            using type = typename base_traits<E1>::type;
            using value_at_t0 = typename base_traits<E1>::value_at_t0;
            using derive_t = typename base_traits<E1>::derive_t;
            using derive_x = typename base_traits<E1>::derive_x;
        };

        // a+x = x+a
        // x not constant
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                    std::enable_if_t<
                        internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value &&
                        !internal::is_instantiation_of<ConstantExpression, typename base_traits<E2>::type>::value 
        >> : base_traits<AddExpression<typename base_traits<E2>::type, typename base_traits<E1>::type>> {};

        // a+b = (a+b)
        // two constants
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                    std::enable_if_t<
                        internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value &&
                        internal::is_instantiation_of<ConstantExpression, typename base_traits<E2>::type>::value
        >> : base_traits<ConstantExpression<
                Q64::add_t<typename base_traits<E1>::type::v, typename base_traits<E2>::type::v>
            >> {};

        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>> {
            using type = SubExpression<E1, E2>;
            using value_at_t0 = FPQ64X::template sub_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = SubExpression<typename base_traits<E1>::derive_t, typename base_traits<E2>::derive_t>;
            using derive_x = SubExpression<typename base_traits<E1>::derive_x, typename base_traits<E2>::derive_x>;
        };

        // x *  y
        // x != 0
        // x != 1
        // y not constant
        // y != (z + w)
        // x != (x+w)
        // y != z-w
        template<typename E1, typename E2>
        struct base_traits<
                MulExpression<E1, E2>,
                std::enable_if_t<
                    !std::is_same<ConstantExpression<typename Q64::zero>, typename base_traits<E1>::type>::value &&
                    !std::is_same<ConstantExpression<typename Q64::one>, typename base_traits<E1>::type>::value &&
                    !internal::is_instantiation_of<ConstantExpression, typename base_traits<E2>::type>::value &&
                    !internal::is_instantiation_of<AddExpression, typename base_traits<E2>::type>::value &&
                    !internal::is_instantiation_of<SubExpression, typename base_traits<E2>::type>::value &&
                    !internal::is_instantiation_of<AddExpression, typename base_traits<E1>::type>::value 
                >
            > {
            using type = MulExpression<E1, E2>;
            using value_at_t0 = FPQ64X::template mul_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = AddExpression<
                            MulExpression<typename base_traits<E1>::derive_t, E2>,
                            MulExpression<E1, typename base_traits<E2>::derive_t>>;
            
            using derive_x = AddExpression<
                            MulExpression<typename base_traits<E1>::derive_x, E2>,
                            MulExpression<E1, typename base_traits<E2>::derive_x>>;
        };

        // 0 * x -> 0
        template<typename E1, typename E2>
        struct base_traits<
                MulExpression<E1, E2>,
                std::enable_if_t<std::is_same<ConstantExpression<typename Q64::zero>,
                                               typename base_traits<E1>::type>::value>> {
            using type = ConstantExpression<typename Q64::zero>;
            using value_at_t0 = typename FPQ64X::zero;
            using derive_t = ConstantExpression<typename Q64::zero>;
            using derive_x = ConstantExpression<typename Q64::zero>;
        };

        // 1 * x -> x
        template<typename E1, typename E2>
        struct base_traits<
                MulExpression<E1, E2>,
                std::enable_if_t<std::is_same<ConstantExpression<typename Q64::one>,
                                               typename base_traits<E1>::type>::value>> {
            using type = typename base_traits<E2>::type;
            using value_at_t0 = typename base_traits<E2>::value_at_t0;
            using derive_t = typename base_traits<E2>::derive_t;
            using derive_x = typename base_traits<E2>::derive_t;
        };

        // a*b -> ab
        template<typename E1, typename E2>
        struct base_traits<
                MulExpression<E1, E2>,
                std::enable_if_t<
                        internal::is_instantiation_of<ConstantExpression, typename base_traits<E2>::type>::value &&
                        internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value &&
                        !std::is_same<typename base_traits<E1>::type, ConstantExpression<typename Q64::zero>>::value &&
                        !std::is_same<typename base_traits<E1>::type, ConstantExpression<typename Q64::one>>::value
                >
        > : ConstantExpression<Q64::mul_t<typename base_traits<E1>::type::v, typename base_traits<E2>::type::v>> {};

        // x * a -> a * x
        // a constant
        // x not constant
        // x != a
        template<typename E1, typename E2>
        struct base_traits<
                MulExpression<E1, E2>,
                std::enable_if_t<
                        internal::is_instantiation_of<ConstantExpression, typename base_traits<E2>::type>::value &&
                        !internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value 
        >> : 
        base_traits<MulExpression<typename base_traits<E2>::type, typename base_traits<E1>::type>> {};

         // x * (a * y) -> a * x
        // a constant not zero or 1
        // x not constant
        // y not constant
        // x not addition or sub
        // y not addition or sub
        template<typename E1, typename E2, typename E3>
        struct base_traits<
                MulExpression<E1, MulExpression<E2, E3>>,
                std::enable_if_t<
                        internal::is_instantiation_of<ConstantExpression, typename base_traits<E2>::type>::value &&
                        !std::is_same<ConstantExpression<typename Q64::zero>, typename base_traits<E2>::type>::value &&
                        !std::is_same<ConstantExpression<typename Q64::one>, typename base_traits<E2>::type>::value &&
                        !internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value &&
                        !internal::is_instantiation_of<ConstantExpression, typename base_traits<E3>::type>::value &&
                        !internal::is_instantiation_of<AddExpression, typename base_traits<E1>::type>::value &&
                        !internal::is_instantiation_of<AddExpression, typename base_traits<E3>::type>::value
        >> : 
            base_traits<
                MulExpression<
                    typename base_traits<E2>::type, 
                    MulExpression<
                        typename base_traits<E1>::type, 
                        typename base_traits<E3>::type
                    >
                >
            > {};

        // x * (y + z)
        // y+z not constant
        // x not constant
        template<typename E1, typename E2, typename E3>
        struct base_traits<
                MulExpression<E1, AddExpression<E2, E3>>, 
                std::enable_if_t<
                    !internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value && 
                    !internal::is_instantiation_of<ConstantExpression, typename base_traits<AddExpression<E2, E3>>::type>::value
                >> : 
            base_traits<
                AddExpression<
                    MulExpression<typename base_traits<E1>::type, typename base_traits<E2>::type>,
                    MulExpression<typename base_traits<E1>::type, typename base_traits<E3>::type>
                >
            > {};

        // (x + y) * z
        // x+y not constant
        // z not constant
        // z not an addition
        template<typename E1, typename E2, typename E3>
        struct base_traits<
                MulExpression<AddExpression<E1, E2>, E3>, 
                std::enable_if_t<
                    !internal::is_instantiation_of<ConstantExpression, typename base_traits<E3>::type>::value && 
                    !internal::is_instantiation_of<AddExpression, typename base_traits<E3>::type>::value && 
                    !internal::is_instantiation_of<ConstantExpression, typename base_traits<AddExpression<E1, E2>>::type>::value
                >> : 
            base_traits<
                AddExpression<
                    MulExpression<typename base_traits<E1>::type, typename base_traits<E3>::type>,
                    MulExpression<typename base_traits<E2>::type, typename base_traits<E3>::type>
                >
            > {};

        // x * (y - z)
        // y-z not constant
        // 
        // x not an addition
        // x not an sub
        template<typename E1, typename E2, typename E3>
        struct base_traits<
                MulExpression<E1, SubExpression<E2, E3>>, 
                std::enable_if_t<
                    //!internal::is_instantiation_of<ConstantExpression, typename base_traits<E1>::type>::value && 
                    !internal::is_instantiation_of<AddExpression, typename base_traits<E1>::type>::value && 
                    //!internal::is_instantiation_of<SubExpression, typename base_traits<E1>::type>::value && 
                    !internal::is_instantiation_of<ConstantExpression, typename base_traits<SubExpression<E2, E3>>::type>::value
                >> : 
            base_traits<
                SubExpression<
                    MulExpression<typename base_traits<E1>::type, typename base_traits<E2>::type>,
                    MulExpression<typename base_traits<E1>::type, typename base_traits<E3>::type>
                >
            > {};

        // // (x - y) * z
        // // x-y not constant
        // // z not constant
        // // z not an substraction
        // // z not an addition
        // template<typename E1, typename E2, typename E3>
        // struct base_traits<
        //         MulExpression<SubExpression<E1, E2>, E3>, 
        //         std::enable_if_t<
        //             !internal::is_instantiation_of<ConstantExpression, typename base_traits<E3>::type>::value && 
        //             !internal::is_instantiation_of<SubExpression, typename base_traits<E3>::type>::value && 
        //             !internal::is_instantiation_of<AddExpression, typename base_traits<E3>::type>::value && 
        //             !internal::is_instantiation_of<ConstantExpression, typename base_traits<SubExpression<E1, E2>>::type>::value
        //         >> : 
        //     base_traits<
        //         SubExpression<
        //             MulExpression<typename base_traits<E1>::type, typename base_traits<E3>::type>,
        //             MulExpression<typename base_traits<E2>::type, typename base_traits<E3>::type>
        //         >
        //     > {};

        template<typename E1, typename E2>
        struct base_traits<DivExpression<E1, E2>> {
            using type = DivExpression<E1, E2>;
            using value_at_t0 = typename FPQ64X::template div_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = DivExpression<
                                SubExpression<
                                    MulExpression<typename base_traits<E1>::derive_t, E2>,
                                    MulExpression<E1, typename base_traits<E2>::derive_t>
                                >,
                                MulExpression<E2, E2>
                            >;
            
            using derive_x = DivExpression<
                                SubExpression<
                                    MulExpression<typename base_traits<E1>::derive_x, E2>,
                                    MulExpression<E1, typename base_traits<E2>::derive_x>
                                >,
                                MulExpression<E2, E2>
                            >;
        };

        template<typename E>
        struct base_traits<ExpExpression<E>> {
            using type = ExpExpression<E>;
            using value_at_t0 = typename FPQ64X::one; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<typename base_traits<E>::derive_t, ExpExpression<E>>;
            using derive_x = MulExpression<typename base_traits<E>::derive_x, ExpExpression<E>>;
        };
        
        template<typename E>
        struct base_traits<SinExpression<E>> {
            using type = SinExpression<E>;
            using value_at_t0 = typename FPQ64X::zero; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<typename base_traits<E>::derive_t, CosExpression<E>>;
            using derive_x = MulExpression<typename base_traits<E>::derive_x, CosExpression<E>>;
        };
        
        template<typename E>
        struct base_traits<CosExpression<E>> {
            using type = CosExpression<E>;
            using value_at_t0 = typename FPQ64X::one; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<typename base_traits<E>::derive_t, SubExpression<ConstantExpression<Q64::zero>, SinExpression<E>>>;
            using derive_x = MulExpression<typename base_traits<E>::derive_x, SubExpression<ConstantExpression<Q64::zero>, SinExpression<E>>>;
        };

        template<typename E>
        struct base_traits<ExpExpression<E>, ExpExpression<SubExpression<ConstantExpression<Q64::zero>, E>>> {
            using type = ConstantExpression<typename Q64::one>;
            using value_at_t0 = typename FPQ64X::one;
            using derive_t = ConstantExpression<typename Q64::zero>;
            using derive_x = ConstantExpression<typename Q64::zero>;
        };

        template<typename E, size_t k>
        struct taylor_coeff_helper {
            using type = FPQ64X::remove_divisor_t<FPQ64X::template div_t<typename taylor_coeff_helper<typename base_traits<E>::type::derive_t, k-1>::type,
                                                FPQ64X::inject_constant_t<k>>>;
        };

        template<typename E>
        struct taylor_coeff_helper<E, 0> {
            using type = typename base_traits<E>::type::value_at_t0;
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

        using T = aad::TExpression<Q64::one, 1>;
        using X = aad::XExpression<Q64::one, 1>;
        using XT = aad::MulExpression<X, T>;
        using ONE = aad::ConstantExpression<Q64::one>;
    }
}