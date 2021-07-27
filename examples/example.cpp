#include "qsl/qubits.hpp"
#include "qsl/utils/random.hpp"
#include "qsl/utils/quantum.hpp"

int main ()
{
    unsigned nqubits = 4;

    std::vector<qsl::complex<double>> rand = qsl::makeRandomNPState<double>(nqubits);
    qsl::Qubits<qsl::Type::NP, double> q{rand};
    qsl::Qubits<qsl::Type::OmpNP, double> qomp{rand};
    qsl::Qubits<qsl::Type::Default, double> qref{rand};
    
    q.pauliZ(0);
    q.print();

    qomp.pauliZ(0);
    qomp.print();

    qref.pauliZ(0);
    qref.print();

    std::cout << "q and qomp = "
	      << qsl::fubiniStudy(q.getState(), qomp.getState())
	      << std::endl;

    
    std::cout << "q and qref = "
	      << qsl::fubiniStudy(q.getState(), qref.getState())
	      << std::endl;

}
