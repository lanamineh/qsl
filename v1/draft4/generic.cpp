#include <iostream>

template<typename... Args>
struct PrintFn
{
    /**
     * \brief Print the types in Args...
     */
    static constexpr void print()
	{
	    ((std::cout << typeid(Args).name() << ","), ...);
	    std::cout << std::endl;
	}
};

template<typename... Args> struct Typelist;

template<>
struct Typelist<> : PrintFn<>
{
    static constexpr std::size_t size = 0;
};

template<typename First, typename... Rest>
struct Typelist<First, Rest...> : PrintFn<First, Rest...>
{
    static constexpr std::size_t size = Typelist<Rest...>::size + 1;
};

int main()
{
    std::cout << Typelist<int,double,void>::size << std::endl;
    Typelist<int,double,void*>::print();
}
