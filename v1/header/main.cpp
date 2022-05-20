#include "qsl.hpp"
#include <iostream>

int main()
{
    qsl::basic<double, true> q{2};
    std::cout << q[10] << std::endl;
}
