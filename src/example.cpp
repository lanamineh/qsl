#include "qsl/qubits.hpp"
#include "qsl/utils/random.hpp"
#include "qsl/utils/quantum.hpp"

int main ()
{
    unsigned nqubits = 4;

    std::vector<qsl::complex<double>> rand = qsl::makeRandomState<double>(nqubits);
    qsl::Qubits<qsl::Type::Default, double> q{rand};
    qsl::Qubits<qsl::Type::Omp, double> qomp{rand};
    
    q.rotateX(0,12);
    q.print();

    qomp.rotateX(0,12);
    qomp.print();

}
