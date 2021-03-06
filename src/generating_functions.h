namespace aerobus {
    // generating function - first specy
    struct Tchebychev1 {
        private:
            using TWO = aad::ConstantExpression<Q64::inject_constant_t<2>>;
            using B0 = aad::DivExpression<
                // 1 - xt
                aad::SubExpression<aad::ONE, aad::XT>,
                // 1 - 2xt + t^2
                aad::AddExpression<type_list<
                    // 1 - 2xt
                    aad::SubExpression<aad::ONE, aad::MulExpression<type_list<TWO, aad::XT>>>,
                    aad::TExpression<Q64::one, 2>>
                >
            >;

        public:
            template<size_t k>
            using at = aad::taylor_coeff_t<B0, k>;
    };

    struct Tchebychev2 {
        private:
        using TWO = aad::ConstantExpression<Q64::inject_constant_t<2>>;
        using B0 = aad::DivExpression<
            // 1 
            aad::ONE,
            // 1 - 2xt + t^2
            aad::AddExpression<type_list<
                aad::SubExpression<aad::ONE, aad::MulExpression<type_list<TWO, aad::XT>>>,
                aad::TExpression<Q64::one, 2>
            >>
        >;

        public:
            template<size_t k>
            using at = aad::taylor_coeff_t<B0, k>;
    };

    struct Euler {
        private:
            using TWO = aad::ConstantExpression<Q64::inject_constant_t<2>>;
            using B0 = aad::DivExpression<
                aad::MulExpression<type_list<
                    TWO,
                    aad::ExpExpression<aad::XT>
                >>,
                aad::AddExpression<type_list<
                    aad::ExpExpression<aad::T>,
                    aad::ONE
                >>
            >;
        using P = polynomial<Q64, 'x'>;
        public:
            template<size_t k>
            using at = P::mul_t<
                            aad::taylor_coeff_t<B0, k>, 
                            typename P::template val<
                                typename Q64::template val<
                                    typename factorial<i64, k>::type, 
                                    typename i64::one
                                >
                            >
                        >;
    };

    template<size_t m>
    struct Bernstein {
    private:
        using P = polynomial<Q64, 'x'>;
        using OneMU = typename P::template val<Q64::inject_constant_t<-1>, typename Q64::one>;
        template<size_t i, size_t index>
        struct inner {
            using type = P::mul_t<OneMU, typename inner<i, index-1>::type>;
        };

        template<size_t i>
        struct inner<i, 0> {
            using type = P::mul_t<
                                P::monomial_t<typename Q64::one, i>, 
                                typename P::template val<
                                    typename Q64::template val<
                                        typename combination<i64, i, m>::type,
                                        typename i64::one
                                    >
                                >
                            >;
        };
    public:

        template<size_t k>
        using at = typename inner<k, m-k>::type;
    };

    template<size_t k>
    struct Hermite {
        private:
        using B0 = aad::ExpExpression<
                aad::SubExpression<
                    aad::ConstantExpression<typename Q64::zero>,
                    aad::MulExpression<type_list<
                        aad::T,
                        aad::T
                    >>
                >
            >;

        using B1 = aad::ExpExpression<
                    aad::MulExpression<type_list<
                        aad::T,
                        aad::T
                    >>
            >;

        template<size_t m, bool dummy = true>
        struct inner {
            using type = typename inner<m-1, dummy>::type::derive_t;
        };

        template<bool dummy>
        struct inner<0, dummy> {
            using type = B0;
        };

    public:
        using at = typename aad::MulExpression<type_list<typename inner<k>::type, B1>>::type;
    };
}