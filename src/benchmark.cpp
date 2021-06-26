#include "qsl/qubits.hpp"
#include "qsl/benchmark/compare.hpp"

int main ()
{
    unsigned nqubits = 3;
    qsl::Qubits<qsl::Type::Default, double> q{nqubits};

    qsl::Compare<qsl::Test::SingleSim,
		 qsl::Qubits<qsl::Type::Default, double>,
		 qsl::Qubits<qsl::Type::Omp, double>> cmp{16,1000};
    cmp.phase();
    
}
