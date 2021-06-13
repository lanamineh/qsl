#include "qsl/qubits.hpp"
#include "qsl/benchmark/compare.hpp"

int main ()
{
    unsigned nqubits = 3;
    Qubits<Type::Default, double> q{nqubits};

    Compare<Test::SingleSim,
	    Qubits<Type::Default, double>,
	    Qubits<Type::Omp, double>> cmp{16,1000};
    cmp.phase();
    
}
