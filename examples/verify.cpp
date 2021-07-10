#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::NP, double>,
		qsl::NPStateGen<double>,
		qsl::NPGateChecker> verify;
    verify.configureState(4, 2);
    auto result = verify.check<qsl::NPGateChecker>();

    for (const auto & [gate,table] : result) {
	std::cout << "Gate name: " << gate << std::endl;
        table.print();
    }
    
}
