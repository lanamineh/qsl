#include <gtest/gtest.h>
#include <qsl/qubits.hpp>
#include <armadillo>
#include <complex>
#include <list>

template<typename T>
struct SimWrapper
{
    using Sim = T;
};

/**
 * \brief Typed test suite for one-qubit gates
 */
template <typename T>
class OneQubitGates : public testing::Test {
public:
    using List = std::list<T>;
    static T shared_;
    T value_;
};

template <typename T>
class NPOneQubitGates : public testing::Test {
public:
    using List = std::list<T>;
    static T shared_;
    T value_;
};


/// List of simulator types to check
using Sim1 = SimWrapper<qsl::Qubits<qsl::Type::Default, float>>;
using Sim2 = SimWrapper<qsl::Qubits<qsl::Type::Default, double>>;
using Sim3 = SimWrapper<qsl::Qubits<qsl::Type::Omp, float>>;
using Sim4 = SimWrapper<qsl::Qubits<qsl::Type::Omp, double>>;
using Sim5 = SimWrapper<qsl::Qubits<qsl::Type::Resize, float>>;
using Sim6 = SimWrapper<qsl::Qubits<qsl::Type::Resize, double>>;
using SimTypes = ::testing::Types<Sim1, Sim2, Sim3, Sim4, Sim5, Sim6>;

using Sim7 = SimWrapper<qsl::Qubits<qsl::Type::NP, float>>;
using Sim8 = SimWrapper<qsl::Qubits<qsl::Type::NP, double>>;
using NPSimTypes = ::testing::Types<Sim7,Sim8>;

TYPED_TEST_SUITE(OneQubitGates, SimTypes);

TYPED_TEST(OneQubitGates, OneQubitNoArg)
{   
    const unsigned num_qubits = 8;
    const unsigned targ = 4;
    using Sim = TypeParam::Sim;//qsl::Qubits<qsl::Type::Default, Fp>; 
    using Fp = TypeParam::Sim::Fp_type;

   
    // Create list of gates mapped to matrices
    std::vector<std::pair<std::function<void(Sim &, unsigned)>, 
			  arma::SpMat<std::complex<Fp>>>> gates;

    // PauliX 
    arma::SpMat<std::complex<Fp>> pauliX(2, 2);
    pauliX(0, 1) = 1;
    pauliX(1, 0) = 1;
    auto fn_pauliX = [](Sim & sim, unsigned targ) {
			 sim.pauliX(targ);
		     };
    gates.push_back({fn_pauliX, pauliX});

    // PauliZ
    arma::SpMat<std::complex<Fp>> pauliZ(2, 2);
    pauliZ(0, 0) = 1;
    pauliZ(1, 1) = -1;
    auto fn_pauliZ = [](Sim & sim, unsigned targ) {
			 sim.pauliZ(targ);
		     };
    gates.push_back({fn_pauliZ, pauliZ});
    
    // Hadamard
    arma::SpMat<std::complex<Fp>> hadamard(2, 2);
    Fp sqrt2 = 1/std::sqrt(2);
    hadamard(0, 0) = sqrt2;
    hadamard(0, 1) = sqrt2;
    hadamard(1, 0) = sqrt2;
    hadamard(1, 1) = -sqrt2;
    auto fn_hadamard = [](Sim & sim, unsigned targ) {
			   sim.hadamard(targ);
		       };
    gates.push_back({fn_hadamard, hadamard});

    // phase shift
    Fp angle = 0.4;
    arma::SpMat<std::complex<Fp>> phase(2, 2);
    phase(0, 0) = 1;
    phase(1, 1) = std::complex<Fp>{std::cos(angle), std::sin(angle)};
    auto fn_phase = [=](Sim & sim, unsigned targ) {
			sim.phase(targ, angle);
		       };
    gates.push_back({fn_phase, phase});

    // rotateZ    
    arma::SpMat<std::complex<Fp>> rotateZ(2, 2);
    rotateZ(0, 0) = std::complex<Fp>{std::cos(angle/2), std::sin(-angle/2)};
    rotateZ(1, 1) = std::complex<Fp>{std::cos(angle/2), std::sin(angle/2)};
    auto fn_rotateZ = [=](Sim & sim, unsigned targ) {
			sim.rotateZ(targ, angle);
		       };
    gates.push_back({fn_rotateZ, rotateZ});
    
    // rotateX
    arma::SpMat<std::complex<Fp>> rotateX(2, 2);
    rotateX(0, 0) = std::cos(angle/2);
    rotateX(0, 1) = std::complex<Fp>{0,-std::sin(angle/2)};
    rotateX(1, 0) = std::complex<Fp>{0,-std::sin(angle/2)};
    rotateX(1, 1) = std::cos(angle/2);
    auto fn_rotateX = [=](Sim & sim, unsigned targ) {
			sim.rotateX(targ, angle);
		       };
    gates.push_back({fn_rotateX, rotateX});
    
    
    for (const auto & [fn, mat] : gates) {    
	// Make a random state
	Sim q{num_qubits};
	const std::vector<qsl::complex<Fp>> state
	    = qsl::makeRandomState<Fp>(num_qubits);
	q.setState(state);

	// Set an armadillo vector to the same state
	std::size_t dim = 1 << num_qubits;
	arma::Col<std::complex<Fp>> v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    v(i) = std::complex<Fp>{state[i].real, state[i].imag};
	}

	// Apply gate to qubits
	std::invoke(fn, q, targ);

	// Create gate in armadillo
	// Sizes of idenity matrix padding
	std::size_t pre = 1 << targ;
	std::size_t post = 1 << (num_qubits - targ - 1);
	// Tensor to make the gate
	arma::SpMat<std::complex<Fp>> gate =
	    arma::speye<arma::SpMat<std::complex<Fp>>>(pre, pre);
	gate = arma::kron(mat, gate);
	gate = arma::kron(arma::speye<arma::SpMat<std::complex<Fp>>>(post, post),
			  gate);

	// Apply gate in armadillo
	v = gate * v;

	// Read qubit state into armadillo
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}
    
    	EXPECT_TRUE(arma::approx_equal(v, qubit_v, "both", 1e-6, 1e-8));
    }
}
