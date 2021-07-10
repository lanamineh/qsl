#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::Omp, double>,
		qsl::DefaultStateGen<double>,
		qsl::DefaultGateChecker> verify;
    verify.configureState(4);
    auto result = verify.check<qsl::DefaultGateChecker>();

    for (const auto & [gate,results] : result.distances) {
	std::cout << "Gate name: " << gate << std::endl;
	results.print();
    }
    
}
