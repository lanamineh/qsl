#include <iostream>

template<typename... Args> struct Typelist;

template<>
struct Typelist<>
{
    static void print() { std::cout << std::endl; };
};

template<typename First, typename... Rest>
struct Typelist<First, Rest...>
{
    static void print() {
	std::cout << typeid(First).name() << ",";
	Typelist<Rest...>::print();
    }
};

int main()
{
    Typelist<int,double,double,unsigned>::print();
}
