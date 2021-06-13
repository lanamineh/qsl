#include "qsl/qubits.hpp"

int main ()
{
    unsigned nqubits = 3;
    Qubits<Type::Default, double> q{nqubits};
    q.pauliX(2);
    q.print();
}
