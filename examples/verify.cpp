#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::Omp, double>,
		qsl::DefaultStateGen<double>,
		qsl::MeasureChecker> verify;
    verify.configureState(4);
    verify.configureChecker<qsl::MeasureChecker>(10000, 0.95);
    auto result = verify.check<qsl::MeasureChecker>();
    //verify.checkAll();
    result.print();
    
}
