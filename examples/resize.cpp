#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    using Sim = qsl::Qubits<qsl::Type::Resize, double>;    

    // Make a simulator in computational basis state
    Sim q{ 4 };
    q.print();
    
    std::vector<qsl::complex<double>> state = q.getState();

    unsigned outcome = q.measureOut(0);
    std::cout << outcome << std::endl;

    q.print();
       
}
