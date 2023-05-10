#ifndef __DEF_TYPE_LIST__
#define __DEF_TYPE_LIST__

namespace aerobus
{
    template <typename... Ts>
    struct type_list;

    namespace internal
    {
        template <typename T, typename... Us>
        struct pop_front_h
        {
            using tail = type_list<Us...>;
            using head = T;
        };

        template <uint64_t index, typename L1, typename L2>
        struct split_h
        {
        private:
            static_assert(index <= L2::length, "index ouf of bounds");
            using a = typename L2::pop_front::type;
            using b = typename L2::pop_front::tail;
            using c = typename L1::template push_back<a>;

        public:
            using head = typename split_h<index - 1, c, b>::head;
            using tail = typename split_h<index - 1, c, b>::tail;
        };

        template <typename L1, typename L2>
        struct split_h<0, L1, L2>
        {
            using head = L1;
            using tail = L2;
        };

        template <uint64_t index, typename L, typename T>
        struct insert_h
        {
            static_assert(index <= L::length, "index ouf of bounds");
            using s = typename L::template split<index>;
            using left = typename s::head;
            using right = typename s::tail;
            using ll = typename left::template push_back<T>;
            using type = typename ll::template concat<right>;
        };

        template <uint64_t index, typename L>
        struct remove_h
        {
            using s = typename L::template split<index>;
            using left = typename s::head;
            using right = typename s::tail;
            using rr = typename right::pop_front::tail;
            using type = typename left::template concat<rr>;
        };
    }

    template <typename... Ts>
    struct type_list
    {
    private:
        template <typename T>
        struct concat_h;

        template <typename... Us>
        struct concat_h<type_list<Us...>>
        {
            using type = type_list<Ts..., Us...>;
        };

    public:
        static constexpr size_t length = sizeof...(Ts);

        template <typename T>
        using push_front = type_list<T, Ts...>;

        template <uint64_t index>
        using at = internal::type_at_t<index, Ts...>;

        struct pop_front
        {
            using type = typename internal::pop_front_h<Ts...>::head;
            using tail = typename internal::pop_front_h<Ts...>::tail;
        };

        template <typename T>
        using push_back = type_list<Ts..., T>;

        template <typename U>
        using concat = typename concat_h<U>::type;

        template <uint64_t index>
        struct split
        {
        private:
            using inner = internal::split_h<index, type_list<>, type_list<Ts...>>;

        public:
            using head = typename inner::head;
            using tail = typename inner::tail;
        };

        template <uint64_t index, typename T>
        using insert = typename internal::insert_h<index, type_list<Ts...>, T>::type;

        template <uint64_t index>
        using remove = typename internal::remove_h<index, type_list<Ts...>>::type;
    };

    template <>
    struct type_list<>
    {
        static constexpr size_t length = 0;

        template <typename T>
        using push_front = type_list<T>;

        template <typename T>
        using push_back = type_list<T>;

        template <typename U>
        using concat = U;

        // TODO: assert index == 0
        template <uint64_t index, typename T>
        using insert = type_list<T>;
    };
}

#endif