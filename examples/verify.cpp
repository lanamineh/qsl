#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::Omp, double>,
		qsl::DefaultStateGen<double>,
		qsl::SampleAllChecker> verify;
    verify.configureState(1);
    verify.configureChecker<qsl::SampleAllChecker>(1000000, 0.99);
    auto result = verify.check<qsl::SampleAllChecker>();
    //verify.checkAll();
    result.print();
    
}
