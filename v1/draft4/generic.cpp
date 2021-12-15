#include <iostream>
#include <tuple>

template<typename... Args>
struct Typelist
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

int main()
{
    using T = Typelist<int,double,void>;

    std::cout << T::size << std::endl;
    T::print();
    
    using R = T::get<1>; 
    std::cout << typeid(R).name() << std::endl;

}
