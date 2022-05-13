#include "sim.hpp"
#include <gtest/gtest.h>

TEST(SimTest, Constructor)
{
    Sim<double> s1{2};
    EXPECT_EQ(s1.nqubits, 2);
    EXPECT_EQ(s1.dim, 4);
    EXPECT_EQ(s1.state.size(), 4);
}

TEST(SimTest, ConstructorException)
{
    EXPECT_THROW(Sim<double> s1{11}, std::logic_error);
}

///////////////////////////////////////////////////////////

class SimTestFix : public ::testing::Test
{
protected:
    Sim<double> s1{2};
};

TEST_F(SimTestFix, Constructor)
{
    EXPECT_EQ(s1.nqubits, 2);
    EXPECT_EQ(s1.dim, 4);
    EXPECT_EQ(s1.state.size(), 4);
}

//////////////////////////////////////////////////////////

template <typename T>
class SimTestType : public testing::Test {};

using MyTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(SimTestType, MyTypes);

TYPED_TEST(SimTestType, Constructor)
{
    using F = TypeParam;

    Sim<F> s1{2};
    EXPECT_EQ(s1.nqubits, 2);
    EXPECT_EQ(s1.dim, 4);
    EXPECT_EQ(s1.state.size(), 4);
}

////////////////////////////////////////////////////

template <std::floating_point F1, bool D1>
struct TypeDefinitions
{
    using F = F1;
    using D = std::bool_constant<D1>;
};

template <typename T>
class SimTestType2 : public testing::Test {};

using MyTypes2 = ::testing::Types<TypeDefinitions<float, true>,
				  TypeDefinitions<float, false>,
				  TypeDefinitions<double, true>,
				  TypeDefinitions<double, false>>;

TYPED_TEST_SUITE(SimTestType2, MyTypes2);

TYPED_TEST(SimTestType2, Constructor)
{
    using F = TypeParam::F;
    constexpr bool D = TypeParam::D::value;
    
    Sim2<F, D> s1{2};
    EXPECT_EQ(s1.nqubits, 2);
    EXPECT_EQ(s1.dim, 4);
    EXPECT_EQ(s1.state.size(), 4);
}

////////////////////////////////////////////////////////////

template <typename T>
struct SimWrapper
{
    using type = T;
};

template <typename T>
class SimTestTypeWrapper : public testing::Test {};

using MyTypesWrapper = ::testing::Types<SimWrapper<Sim2<float, true>>,
					SimWrapper<Sim2<float, false>>>;

TYPED_TEST_SUITE(SimTestTypeWrapper, MyTypesWrapper);

TYPED_TEST(SimTestTypeWrapper, Constructor)
{
    using SimType = TypeParam::type; 
    
    SimType s1{2};
    EXPECT_EQ(s1.nqubits, 2);
    EXPECT_EQ(s1.dim, 4);
    EXPECT_EQ(s1.state.size(), 4);
}


/////////////////////////////////////////////////////////////
template<std::floating_point F, class D>
struct Case
{
    using type = Sim2<F, D::value>;
};

template<class Tuple1, class Tuple2, std::size_t I>
struct make_case
{
    static constexpr std::size_t N = std::tuple_size<Tuple2>::value;

    using type = Case<typename std::tuple_element<I / N, Tuple1>::type,
                      typename std::tuple_element<I % N, Tuple2>::type>;
};

template <class Is, class T1, class... T2>
struct make_combinations;

template <class Tuple1, class Tuple2, std::size_t... Is>
struct make_combinations<std::index_sequence<Is...>, Tuple1, Tuple2>
{
    using tuples = std::tuple<typename make_case<Tuple1, Tuple2, Is>::type...>;
};

template<class Tuple1, class Tuple2>
using Combinations_t = typename make_combinations
                       <std::make_index_sequence<
			    (std::tuple_size<Tuple1>::value)*(std::tuple_size<Tuple2>::value)>,
			Tuple1, Tuple2
			>::tuples;

template <class T>
class MyTestFixture : public ::testing::Test {};

template<class T>
struct Test;

template<class ...T>
struct Test<std::tuple<T...>>
{
    using Types = ::testing::Types<T...>;
};

using TestTypes = Test<Combinations_t<std::tuple<float, double, long double>,
				      std::tuple<std::true_type, std::false_type>>>::Types;

TYPED_TEST_SUITE(MyTestFixture, TestTypes);

TYPED_TEST(MyTestFixture, Test12)
{
    // X or Y
    // TypeParam::type t;

    // // "a" or "b"
    // const std::string str = TypeParam::GetParam();
    using SimType = TypeParam::type;
    
    SimType s1{2};
    EXPECT_EQ(s1.nqubits, 2);
    EXPECT_EQ(s1.dim, 4);
    EXPECT_EQ(s1.state.size(), 4);

}


//////////////////////////////////////////////////////////////////
template <typename T>
class MyTestFixture2 : public testing::Test {};

// Tuple concatenation
template <typename T, typename...>
struct cat
{
    using type = T;
};

template <typename... T1, typename... T2, typename... T3>
struct cat<std::tuple<T1...>, std::tuple<T2...>, T3...>
{
    using type = cat<std::tuple<T1..., T2...>, T3...>::type;
};

template <typename... T>
using cat_t = cat<T...>::type;


// Prepend a type to the beginning of every tuple in a list of tuples
template <typename T, typename V>
struct prepend{};

// Takes std::tuple<std::tuples> as input
template <typename T, typename... L>
struct prepend<T, std::tuple<L...>>
{
    using type = std::tuple<cat_t<std::tuple<T>, L>...>;
};

template <typename T, typename V>
using prepend_t = prepend<T, V>::type;


// Combinations
template <typename... T>
struct comb {};

template <typename... T>
struct comb<std::tuple<T...>>
{
    using type = std::tuple<std::tuple<T>...>;
};

template <typename... T1, typename... T2>
struct comb<std::tuple<T1...>, T2...>
{
    using combined = comb<T2...>::type;

    using type = cat_t<prepend_t<T1, combined>...>;
};

template <typename... T>
using comb_t = comb<T...>::type;



using TestTypes2 = ::testing::Types<cat_t<std::tuple<double>,
					  std::tuple<float, long double>,
					  std::tuple<int>>>;

//using TestTypes3 = ::testing::Types<head_t<std::tuple<double, float, long double>>>;

using TestTypes4 = Test<comb_t<std::tuple<double, char>,
			       std::tuple<float, long double>,
			       std::tuple<int, unsigned>>>::Types;

//using TestTypes4 = ::testing::Types<comb_t<std::tuple<float, int>, std::tuple<double, char>>>;


using TestTypes5 = ::testing::Types<append_t<double, std::tuple<std::tuple<float>, std::tuple<int>>>>;


TYPED_TEST_SUITE(MyTestFixture2, TestTypes4);

TYPED_TEST(MyTestFixture2, Test12)
{
    EXPECT_EQ(0, 0);
}


