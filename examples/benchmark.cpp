#include "qsl/qubits.hpp"
#include "qsl/benchmark/compare.hpp"
#include "qsl/benchmark/time.hpp"
#include "qsl/verify.hpp"

int main ()
{
    unsigned nqubits = 3;
    qsl::Qubits<qsl::Type::Default, double> q{nqubits};

    using Sim1 = qsl::Qubits<qsl::Type::Experimental, double>;
    using Sim2 = qsl::Qubits<qsl::Type::Default, double>;
    
    qsl::Verify<Sim1, Sim2, qsl::DefaultStateGen<double>,
		qsl::DefaultGateChecker> verify;
    verify.configureState(12);
    auto result = verify.check<qsl::DefaultGateChecker>();
    result.at("phase").print();
	
    qsl::Compare<qsl::Test::SingleSim,Sim1, Sim2> cmp{8,100};
    cmp.phase();
    cmp.pauliX();
    cmp.rotateX();
    cmp.controlNot();

    // using Sim1 = qsl::Qubits<qsl::Type::Default, double>;
   
    // qsl::Time<qsl::Test::MultiSim, Sim1, qsl::Restrictions::None> time(12, 200, 100);
    // qsl::Results res1 = time.pauliX();
    // //res1.writeToFile("pauliX.csv");
    // qsl::Results res2 = time.rotateX();
    // //res2.writeToFile("rotateX.csv");
    // qsl::Results res3 = time.controlNot();
    // //res3.writeToFile("cNot.csv");
    // qsl::Results res4 = time.controlPhase();
    // //res4.writeToFile("cPhaseShift.csv");

    
}
