#include <iostream>
#include "sim.hpp"

int main()
{
    Sim<double> s{2};
    std::cout << s.nqubits << std::endl;
    std::cout << s.dim << std::endl;
    std::cout << s.state.size() << std::endl;
    std::cout << s.state[0] << std::endl;
    
    return 0;
}
