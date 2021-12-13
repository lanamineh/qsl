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

int main()
{
    C<double, A, B> c;
    //std::cout << c.a << std::endl;
    //std::cout << c.b << std::endl;
}
