#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::Omp, double>,
		qsl::DefaultStateGen<double>,
		qsl::SampleChecker> verify;
    verify.configureState(4);
    verify.configureChecker<qsl::SampleChecker>(10000, 0.99);
    auto result = verify.check<qsl::SampleChecker>();
    //verify.checkAll();
    result.print();
    
}
