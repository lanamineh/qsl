#include "typelist.hpp"

struct A_def : Replaces<A_def>
{
    A_def() { std::cout << "I am an A_def" << std::endl; }
};

struct B_def : Replaces<B_def>
{
    B_def() { std::cout << "I am a B_def" << std::endl; }
};

struct C_def : Replaces<C_def>
{
    C_def() { std::cout << "I am a C_def " << std::endl; }
};

struct A : Replaces<A_def>
{
    A() { std::cout << "I am an A" << std::endl; }
};

struct B : Replaces<B_def>
{
    B() { std::cout << "I am a B" << std::endl; }
};

struct C : Replaces<C_def>
{
    C() { std::cout << "I am a C" << std::endl; }
};

template<typename... Args>
struct TypeInheriter : Args...
{
};

template<typename TL>
struct Typelist2TypeInheriter;

template<typename... Args>
struct Typelist2TypeInheriter<Typelist<Args...>>
{
    using type = TypeInheriter<Args...>;
};


template<typename... Args>
struct ParseOptions
{
    using def = Typelist<A_def, B_def, C_def>;
    using opt = Typelist<Args...>;
    using parsed = typename opt::parse<def>;
    
    using type = Typelist2TypeInheriter<parsed>::type;
};



template<typename... Args>
struct Generic : ParseOptions<Args...>::type
{
    Generic() { std::cout << "I am Generic " << std::endl; }    
};


template<typename Av, typename Bv, typename Cv>
struct Generic : Av, Bv, Cv
{
    Generic() { std::cout << "I am Generic " << std::endl; }    
};


template<typename TL>
struct GenericWrapper
{
    using type = Generic<Args...>;    
};









int main()
{
    Generic<B> g;   
}



