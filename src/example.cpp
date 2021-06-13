#include "qsl/qubits.hpp"
#include "qsl/utils/random.hpp"
#include "qsl/utils/quantum.hpp"

int main ()
{
    unsigned nqubits = 4;

    std::vector<complex<double>> rand = makeRandomState<double>(nqubits);
    Qubits<Type::Default, double> q{rand};
    Qubits<Type::Omp, double> qomp{rand};
    
    q.pauliX(2);
    q.print();

    qomp.pauliX(2);
    qomp.print();

}
