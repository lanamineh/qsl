#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::Omp, double>,
		qsl::DefaultStateGen<double>,
		qsl::ProbChecker> verify;
    verify.configureState(4);
    auto result = verify.check<qsl::ProbChecker>();
    //verify.checkAll();
    result.print();
    
}
