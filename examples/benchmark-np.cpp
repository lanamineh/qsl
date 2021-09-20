#include "qsl/qubits.hpp"
#include "qsl/qubits.hpp"
#include "qsl/benchmark/compare.hpp"
#include "qsl/benchmark/time.hpp"

int main ()
{

    unsigned nqubits = 16;
    std::size_t dim = 1 << nqubits;

    using Sim = qsl::Qubits<qsl::Type::Default, double>;
    using SimOmp = qsl::Qubits<qsl::Type::Omp, double>;
    using SimNP = qsl::Qubits<qsl::Type::NP, double>;
    using SimNPOmp = qsl::Qubits<qsl::Type::OmpNP, double>;
	
    qsl::Results<unsigned, double> res{{"nones", "Default", "Default+OMP",
					"NP", "NP+OMP"}};
    res.addMeta("Testing the phase gate");
    
    Sim q(nqubits);
    SimOmp q_omp(nqubits);
    SimNP qNP(nqubits);
    SimNPOmp qNP_omp(nqubits);
    
    // Test every possible number of ones
    for (unsigned nones = 1; nones < nqubits; nones++) {
        // Generate a random state
	std::vector<qsl::complex<double>> state =
	    qsl::makeRandomNPState(nqubits, nones);

	// Create some random phases
	std::vector<double> phases = qsl::makeRandomPhases(nqubits);
	
	// Set the states
	q.setState(state);
	q_omp.setState(state);
	qNP.setState(state);
	qNP_omp.setState(state);

	std::cout << "Testing " << nones << " ones -------------" << std::endl;

	// Default
	qsl::Timer tmr;
	tmr.start();
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned i = 0; i < 1000; i++) {
		q.phase(n, phases[n]);
	    }
	}
	tmr.stop();
	double t0 = tmr.getElapsed() / (nqubits * 1000);

	// OMP
	tmr.start();
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned i = 0; i < 1000; i++) {
		q_omp.phase(n, phases[n]);
	    }
	}
	tmr.stop();
	double t1 = tmr.getElapsed() / (nqubits * 1000);

	// NP
	tmr.start();
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned i = 0; i < 1000; i++) {
		qNP.phase(n, phases[n]);
	    }
	}
	tmr.stop();
	double t2 = tmr.getElapsed() / (nqubits * 1000);

	// NP OMP
	tmr.start();
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned i = 0; i < 1000; i++) {
		qNP_omp.phase(n, phases[n]);
	    }
	}
	tmr.stop();
	double t3 = tmr.getElapsed() / (nqubits * 1000);
	
	res.addRow({nones}, {t0, t1, t2, t3});
	  
    }

    res.writeToFile("np_phase.csv");

}
