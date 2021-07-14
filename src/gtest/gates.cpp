#include <gtest/gtest.h>
#include <qsl/qubits.hpp>
#include <qsl/concepts.hpp>
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

/**
 * \brief A class for holding the 
 *
 *
 */
template<qsl::Simulator Sim, typename... Args>
class GateAndMatrix
{
    /**
     * \brief The gate member function to execute
     *
     * Args is either one or two unsigned arguments,
     * depending if the gate is a one-qubit gate
     * or a two-qubit gate.
     *
     */
    std::function<void(Sim &, Args... args)> gate;

    /**
     * \brief The armadillo matrix for the gate
     * 
     * This member holds the armadillo matrix for the one-
     * or two-qubit gate.
     */
    arma::SpMat<std::complex<typename Sim::Fp_type>> mat;

public:

    /**
     * \brief Return the small armadillo matrix
     */
    arma::SpMat<std::complex<typename Sim::Fp_type>> getMat() const
	{
	    return mat;
	}
    
    using Fp = Sim::Fp_type;
    
    GateAndMatrix(std::function<void(Sim &, Args... args)> gate_in,
		  arma::SpMat<std::complex<Fp>> mat_in)
	: gate{gate_in}, mat{mat_in}
	{}

    /**
     * \brief Run the gate on the simulator
     * 
     * Pass either one or two unsigned arguments, which correspond to
     * the qubit indices of the gate which is applied.
     *
     * \return The vector  
     * 
     */
     std::vector<qsl::complex<Fp>> runGate(Sim & sim, Args... args)
	{
	    std::invoke(sim, gate, args...);
	    return sim.getState();
	}

    /**
     * \brief Run the armadillo matrix multiplication
     *
     * Multiply the big armadillo matrix by the state vector. This 
     * function has a different implementation for one- and two-qubit
     * gates, because in the case of two-qubit gates the
     * state vector must be rearranged before the multiplication.
     */
    virtual arma::Col<std::complex<Fp>>
    runArma(arma::Col<std::complex<Fp>>,
	    unsigned num_qubits, Args... args) = 0;
};


template<qsl::Simulator Sim>
class OneQubit : public GateAndMatrix<Sim,unsigned>
{
    using Fp = GateAndMatrix<Sim,unsigned>::Fp;
    
    /**
     * \brief Create armadillo gate
     *
     * This function creates the big armadillo matrix which
     * gets multipled my the state vector. For the single
     * qubit gate you Tensor product the small gate together
     * 
     */
    arma::SpMat<std::complex<Fp>> makeArmaGate(unsigned targ,
					       unsigned num_qubits) const
	{
	    // Create gate in armadillo
	    // Sizes of idenity matrix padding
	    std::size_t pre = 1 << targ;
	    std::size_t post = 1 << (num_qubits - targ - 1);

	    arma::SpMat<std::complex<typename Sim::Fp_type>> mat{
		this->getMat()
	    };
	    
	    // Tensor to make the gate
	    arma::SpMat<std::complex<Fp>> gate =
		arma::speye<arma::SpMat<std::complex<Fp>>>(pre, pre);
	    gate = arma::kron(mat, gate);
	    gate =
		arma::kron(arma::speye<arma::SpMat<std::complex<Fp>>>(post, post),
			   gate);
	    return gate;
	}
    
public:

    OneQubit(std::function<void(Sim &, unsigned)> gate_in,
	     arma::SpMat<std::complex<Fp>> mat_in)
	: GateAndMatrix<Sim,unsigned>{gate_in, mat_in}
	{}
    
    /**
     * \brief Run the armadillo matrix multiplication
     *
     * Multiply the big armadillo matrix by the state vector.
     */
    arma::Col<std::complex<Fp>> runArma(const arma::Col<std::complex<Fp>> & v,
					unsigned num_qubits, unsigned targ)
	{
	    return makeArmaGate(targ, num_qubits) * v;
	}
};

template<qsl::Simulator Sim>
class PauliX : public OneQubit<Sim>
{

    arma::SpMat<std::complex<Fp>> makeMat() const
	{
	    arma::SpMat<std::complex<Fp>> pauliX(2, 2);
	    pauliX(0, 1) = 1;
	    pauliX(1, 0) = 1;
	    return pauliX;

	}

    auto makeFn() const
	{
	    auto fn = [](Sim & sim, unsigned targ) {
			  sim.pauliX(targ);
		      };
	    return fn;
	}
    
public:
    PauliX()
	: 
    
};

TYPED_TEST(OneQubitGates, OneQubitNoArg)
{   
    const unsigned num_qubits = 8;
    const unsigned targ = 4;
    using Sim = TypeParam::Sim;
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
