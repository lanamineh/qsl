#include <iostream>
#include <tuple>

template<typename... Args>
struct TypelistUtils
{
    /// The length of the type list
    static constexpr std::size_t size = sizeof...(Args);
    
    /**
     * \brief Print the types in Args...
     */
    static constexpr void print()
	{
	    ((std::cout << typeid(Args).name() << ","), ...);
	    std::cout << std::endl;
	}

    /**
     * \brief Get the nth type in the list
     */
    template<std::size_t n>
    using get = typename std::tuple_element<n, std::tuple<Args...>>::type;
};

template<typename... Args> struct Typelist;

template<>
struct Typelist<> : TypelistUtils<>
{
    
};

template<typename First, typename... Rest>
struct Typelist<First, Rest...> : TypelistUtils<First, Rest...>
{
    
};

template<typename TL1, typename TL2>
struct Concat;

template<typename... Args1, typename... Args2>
struct Concat<Typelist<Args1...>, Typelist<Args2...>>
{
    using type = Typelist<Args1..., Args2...>;
};

template<typename TL1, typename TL2>
using Concat_t = Concat<TL1, TL2>::type; 

template<typename TL, typename Find, typename Rep>
struct Replace
{};

template<typename TL, typename Find, typename Rep>
using Replace_t = Replace<TL, Find, Rep>::type;

template<typename Find, typename Rep>
struct Replace<Typelist<>, Find, Rep>
{
    using type = Typelist<>;
};

template<typename F, typename R, typename First, typename... Rest>
struct Replace<Typelist<First, Rest...>, F, R>
{
    static constexpr bool found = std::is_same_v<First, F>;

    using next = Replace_t<Typelist<Rest...>, F, R>;
    using replaced = Concat_t<Typelist<R>, next>;
    using kept = Concat_t<Typelist<First>, next>;
    
    using type = std::conditional_t<found, replaced, kept>;
};


template<typename Default, typename Options>
struct Parse;

template<typename Default, typename Options>
using Parse_t = Parse<Default, Options>::type;

template<typename... ArgsDef>
struct Parse<Typelist<ArgsDef...>, Typelist<> >
{
    using type = Typelist<ArgsDef...>;
};

template<typename Def>
struct Replaces
{
    using def = Def;
};

template<typename T>
concept ContainsDef = requires
{
    typename T::def;
};

template<typename... ArgsDef, ContainsDef First, typename... Rest>
struct Parse<Typelist<ArgsDef...>, Typelist<First, Rest...> >
{
    using res = Replace_t<Typelist<ArgsDef...>, typename First::def, First>; 
    using type = Parse_t<res, Typelist<Rest...>>; 
    //using type = Parse_t<Typelist<ArgsDef...>::replace<typename First::def,
    //                                                   First>,
    //                     Typelist<Rest...>>; 

};

struct A_def : Replaces<A_def>
{
};

struct B_def : Replaces<B_def>
{
};

struct C_def : Replaces<C_def>
{
};

struct A : Replaces<A_def>
{
};

struct B : Replaces<B_def>
{
};

struct C : Replaces<C_def>
{
};


int main()
{
    using T = Typelist<A_def, B_def, C_def>;

    using R = Typelist<A, int>;

    using S = Parse_t<T, R>;

    //using S = Concat<T, R>::type;
    //S::print();
    
    //using R = Replace<Typelist<int, double, unsigned>, unsigned, A>::type;
    // using R = Replace<T, A_def, A>::type;
    
    S::print();
}
