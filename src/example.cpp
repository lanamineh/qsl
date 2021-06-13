#include "qsl/qubits.hpp"
#include "qsl/utils/random.hpp"
#include "qsl/utils/quantum.hpp"

int main ()
{
    unsigned nqubits = 4;

    std::vector<complex<double>> rand = makeRandomState<double>(nqubits);
    Qubits<Type::Default, double> q{rand};
    Qubits<Type::Omp, double> qomp{rand};
    
    q.rotateX(0,12);
    q.print();

    qomp.rotateX(0,12);
    qomp.print();

}
