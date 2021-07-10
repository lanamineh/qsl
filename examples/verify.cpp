#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    qsl::Verify<qsl::Qubits<qsl::Type::Default, double>,
		qsl::Qubits<qsl::Type::NP, double>,
		qsl::NPStateGen<double>,
		qsl::PostselectChecker> verify;
    verify.configureState(4,2);
    //verify.configureChecker<qsl::PostselectChecker>(1000000, 0.99);
    auto result = verify.check<qsl::PostselectChecker>();
    //verify.checkAll();
    result.print();
    
}
