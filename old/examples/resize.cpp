#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"
#include "qsl/utils.hpp"

int main ()
{
    using Sim = qsl::Qubits<qsl::Type::Resize, double>;    

    // Make a simulator in computational basis state
    Sim q{ 4 };

    std::vector<qsl::complex<double>> in_state
	= qsl::makeRandomState(4);
    q.setState(in_state);
    
    q.print();


    
    std::vector<qsl::complex<double>> state = q.getState();

    unsigned outcome = q.measureOut(2);
    std::cout << outcome << std::endl;
    q.print();


    q.appendQubit();
    q.print();

}
