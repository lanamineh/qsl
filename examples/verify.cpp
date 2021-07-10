#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::NP, double>,
		qsl::NPStateGen<double>,
		qsl::NPGateChecker> verify;
    verify.configureState(4,2);
    //verify.configureChecker<qsl::PostselectChecker>(1000000, 0.99);
    auto result = verify.check<qsl::NPGateChecker>();
    //verify.checkAll();
    for (const auto & [gate,table] : result) {
	std::cout << "Gate: " << gate << std::endl;
	table.print();
    }
    
}
