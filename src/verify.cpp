#include "qsl/qubits.hpp"
#include "qsl/verify.hpp"

int main ()
{
    Verify<Qubits<Type::Default, double>,
	   Qubits<Type::Omp, double>,
	   DefaultStateGen<double>,
	   DefaultGateChecker> verify;
    verify.configureState(4);
    verify.checkAll();
    
}
