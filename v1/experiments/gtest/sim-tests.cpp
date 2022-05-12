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

template <class T1, class T2, class Is>
struct make_combinations;

template <class Tuple1, class Tuple2, std::size_t... Is>
struct make_combinations<Tuple1, Tuple2, std::index_sequence<Is...>>
{
    using tuples = std::tuple<typename make_case<Tuple1, Tuple2, Is>::type...>;
};

template<class Tuple1, class Tuple2>
using Combinations_t = typename make_combinations
                       <Tuple1, Tuple2,
			std::make_index_sequence<
			    (std::tuple_size<Tuple1>::value)*(std::tuple_size<Tuple2>::value)>
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
