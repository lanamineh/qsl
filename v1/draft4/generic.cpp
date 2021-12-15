#include <iostream>

template<typename... Args>
struct Typelist
{
    /**
     * \brief Print the types
     *
     */
    static void print()
	{
	    ((std::cout << typeid(Args).name() << ","), ...);
	    std::cout << std::endl;
	}
};

int main()
{
    Typelist<>::print();

}
