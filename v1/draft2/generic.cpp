#include <vector>
#include <iostream>


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

template<typename... Args>
class C : private Args...
{
    
};

int main()
{
    C<A,B> c;
    std::cout << c.a << std::endl;
    std::cout << c.b << std::endl;
}
