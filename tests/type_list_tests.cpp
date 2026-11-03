#include <gtest/gtest.h>
#include <type_traits>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(type_list, basic_assertions) {
    using A = type_list<int, float>;
    EXPECT_EQ(A::length, 2);
    EXPECT_TRUE((std::is_same_v<A::at<0>, int>));
    EXPECT_TRUE((std::is_same_v<A::at<1>, float>));

    using B = A::push_back<double>;
    EXPECT_EQ(B::length, 3);
    EXPECT_TRUE((std::is_same_v<B::at<2>, double>));

    using C = B::pop_front;
    EXPECT_EQ(C::tail::length, 2);
    EXPECT_TRUE((std::is_same_v<C::type, int>));
    EXPECT_TRUE((std::is_same_v<C::tail::at<0>, float>));

    using D = C::tail::split<1>;
    EXPECT_EQ(D::head::length, 1);
    EXPECT_EQ(D::tail::length, 1);
    EXPECT_TRUE((std::is_same_v<D::head::at<0>, float>));
    EXPECT_TRUE((std::is_same_v<D::tail::at<0>, double>));
}

TEST(type_list, more_complex_assertions) {
    {
        using t = type_list<float, int, double>;
        using t2 = t::remove<1>;
        EXPECT_EQ(t2::length, 2);
        EXPECT_TRUE((std::is_same_v<t2::at<0>, float>));
        EXPECT_TRUE((std::is_same_v<t2::at<1>, double>));
    }
    {
        using t = type_list<float, double>;
        using t2 = t::insert<1, int>;
        EXPECT_EQ(t2::length, 3);
        EXPECT_TRUE((std::is_same_v<t2::at<0>, float>));
        EXPECT_TRUE((std::is_same_v<t2::at<1>, int>));
        EXPECT_TRUE((std::is_same_v<t2::at<2>, double>));
    }
    {
        using t = type_list<>;
        using t2 = t::insert<0, float>;
        EXPECT_EQ(t2::length, 1);
        EXPECT_TRUE((std::is_same_v<t2::at<0>, float>));
    }
}
