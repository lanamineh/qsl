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

template<typename TL, typename Find, typename Rep>
struct Replace
{};

template<typename Find, typename Rep>
struct Replace<Typelist<>, Find, Rep>
{
    using type = Typelist<>;
};

template<typename Find, typename Rep, typename First, typename... Rest>
struct Replace<Typelist<First, Rest...>, Find, Rep>
{
    // Compare First with Find.
    // if First == Find, return Typelist<Rep, Rest...>
    // else return Typelist<First, Rest...>

    static constexpr bool found = std::is_same_v<First, Find>;
    
    using type = std::conditional_t<found,
				    Replace<Typelist<Rest...>, Find, Rep>,
				    Typelist<First, Rest...>>;
};





struct A_def{};
struct B_def{};
struct C_def{};

struct A{};
struct B{};
struct C{};


int main()
{
    using T = Typelist<A_def, B_def, C_def>;

    using R = Replace<Typelist<int, double, unsigned>, int, double>::type;
    // using R = Replace<T, A_def, A>::type;
    
    R::print();
}
