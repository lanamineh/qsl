#include <vector>
#include <iostream>
#include <tuple>

class A_def
{
public:
    int a = 3;
    A_def() { std::cout << "A_def" << std::endl; }
};

class B_def
{
public:
    int b = 4;
    B_def() { std::cout << "B_def" << std::endl; }    
};

class A
{
public:
    int a = 3;
    A() { std::cout << "A" << std::endl; }
};

class B
{
public:
    int b = 4;
    B() { std::cout << "B" << std::endl; }    
};

template<std::floating_point Fp, typename... Args>
class C : public Args...
{
    
    
};

template<typename... Others>
struct TypeList
{};

/// General case
template<typename First, typename... Others>
struct TypeList<First, Others...>
{
    static void print()
	{
	    std::cout << typeid(First).name() << ",";
	    TypeList<Others...>::print();
	}
};

/// Base case
template<>
struct TypeList<>
{
    static void print()
	{
	    std::cout << std::endl;
	}
};




template<typename TL_in, typename TL_out = TypeList<>>
struct UniqueTypeList
{};

// Check whether type First is in typenames Others...
template<typename First, typename... Others>
static constexpr bool contains {
    std::disjunction_v<std::is_same<First, Others>...>
};


template<typename First, typename... Others, typename... Pruned>
struct UniqueTypeList<TypeList<First, Others...>, TypeList<Pruned...>>
{    
    // Whether to push the new type of not
    static constexpr bool push = contains<First, Pruned...>;

    // The new type including First
    using WithFirst = TypeList<Pruned..., First>; 

    // The new type excluding First
    using WithoutFirst = TypeList<Pruned...>; 
    
    // The new typelist
    using NewTypeList = std::conditional_t<push, WithoutFirst, WithFirst>;
    
    // Recursive traversal of the list
    using next = UniqueTypeList<TypeList<Others...>, NewTypeList>::next;
};

template<typename... Pruned>
struct UniqueTypeList<TypeList<>, TypeList<Pruned...>>
{
    // The new type excluding First
    using next = TypeList<Pruned...>; 
    
    // No need -- Recursive traversal of the list
    //using next = UniqueTypeList<TypeList<Others...>, NewTypeList>::next;
};

int main()
{
    //C<double, A, B> c;
    //std::tuple<A,B> tuple;
    //std::cout << c.a << std::endl;
    //std::cout << c.b << std::endl;

    using TestTL = TypeList<double,int,float,double,double,int,unsigned,unsigned>;
    
    TestTL::print();

    using UniqueTL = UniqueTypeList<TestTL>::next;
    UniqueTL::print();
}
