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
TYPED_TEST_SUITE(NPOneQubitGates, SimTypes);

/// Gets the bit in the nth position of val
unsigned getBit(unsigned val, unsigned n)
{
    return 1 & (val >> n);
}

/// Set the nth bit of val to b
unsigned setBit(unsigned val, unsigned n, unsigned b)
{
    if (b == 0) {
	val &= ~(1 << n);
    } else {
	val |= (1 << n);
    }
    return val;
}


/// Good luck...
template<std::floating_point Fp>
arma::SpMat<std::complex<Fp>>
makeMatrix(const arma::Mat<std::complex<Fp>> & gate, unsigned nqubits,
	   unsigned ctrl, unsigned targ)
{
    const std::size_t dim{ 1ULL << nqubits };
    arma::SpMat<std::complex<Fp>> mat(dim,dim);    

    // Write columns of mat. The columns of the matrix
    // correspond to the images of the basis states under
    // the action of the unitary matrix. 
    for (std::size_t col = 0; col < dim; col++) {

	// Get column bitstring and find the value of the ctrl and targ bits.
	// The values of the bits of col in the ctrl and targ positions
	// fix which column of the small matrix is used in populating
	// values of the big matrix at column col.
	const unsigned ctrl_val = getBit(col, ctrl);
	const unsigned targ_val = getBit(col, targ);

	// Write rows of mat. For each row index, 
	for (std::size_t n = 0; n < gate.n_rows; n++) {

	    // Make the row index. The rows that are non-zero in a
	    // particular column are the ones whose bits agree with
	    // the bits in the col, apart from at the ctrl and targ
	    // positions. There, they take every possible ctrl and
	    // targ value.
	    std::size_t row = col;
	    row = setBit(row, ctrl, getBit(n,1));
	    row = setBit(row, targ, getBit(n,0));

	    // Make the column index for the small matrix
	    std::size_t k = (ctrl_val << 1) | targ_val;
	    mat(row,col) = gate(n,k);
	}
    }
    return mat;
}

TEST(GateTests,MakeMatrixTestsCnot)
{

    // CNOT gate, qubit 1 is ctrl, qubit 0 is targ
    arma::Mat<std::complex<double>> gate(4, 4, arma::fill::zeros);
    gate(0b00,0b00) = 1;
    gate(0b01,0b01) = 1;
    gate(0b10,0b11) = 1;
    gate(0b11,0b10) = 1;

    // Apply CNOT to qubits ctrl = 1, targ = 0
    unsigned ctrl = 1;
    unsigned targ = 0;
    const auto m0 = makeMatrix(gate, 3, ctrl, targ);
    std::cout << m0 << std::endl;
    // Check all the basis states
    const std::complex<double> one{1,0};
    EXPECT_EQ(m0(0b000,0b000), one); // |000) -> |000) 
    EXPECT_EQ(m0(0b001,0b001), one); // |001) -> |001) 
    EXPECT_EQ(m0(0b010,0b011), one); // |010) -> |011) 
    EXPECT_EQ(m0(0b011,0b010), one); // |011) -> |010) 
    EXPECT_EQ(m0(0b100,0b100), one); // |100) -> |100) 
    EXPECT_EQ(m0(0b101,0b101), one); // |101) -> |101) 
    EXPECT_EQ(m0(0b110,0b111), one); // |110) -> |111) 
    EXPECT_EQ(m0(0b111,0b110), one); // |111) -> |110) 

    // Apply CNOT to qubits ctrl = 0, targ = 2
    ctrl = 0;
    targ = 2;
    const auto m1 = makeMatrix(gate, 3, ctrl, targ);
    std::cout << m1 << std::endl;
    // Check all the basis states
    EXPECT_EQ(m1(0b000,0b000), one); 
    EXPECT_EQ(m1(0b001,0b101), one);
    EXPECT_EQ(m1(0b010,0b010), one);
    EXPECT_EQ(m1(0b011,0b111), one);
    EXPECT_EQ(m1(0b100,0b100), one);
    EXPECT_EQ(m1(0b101,0b001), one);
    EXPECT_EQ(m1(0b110,0b110), one);
    EXPECT_EQ(m1(0b111,0b011), one);

    // Check the degnerate cases
    ctrl = 0;
    targ = 1;
    const auto m2 = makeMatrix(gate, 2, ctrl, targ);
    std::cout << m2 << std::endl;
    // Check all the basis states
    EXPECT_EQ(m2(0b00,0b00), one); 
    EXPECT_EQ(m2(0b01,0b11), one);
    EXPECT_EQ(m2(0b10,0b10), one);
    EXPECT_EQ(m2(0b11,0b01), one);

    // Check the degnerate case the other way round
    ctrl = 1;
    targ = 0;
    const auto m3 = makeMatrix(gate, 2, ctrl, targ);
    std::cout << m3 << std::endl;
    // Check all the basis states
    EXPECT_EQ(m3(0b00,0b00), one); 
    EXPECT_EQ(m3(0b01,0b01), one);
    EXPECT_EQ(m3(0b10,0b11), one);
    EXPECT_EQ(m3(0b11,0b10), one);
}

TEST(GateTests,MakeMatrixTestsSwap)
{

    // SWAP gate, qubit 1 is ctrl, qubit 0 is targ
    arma::Mat<std::complex<double>> gate(4, 4, arma::fill::zeros);
    gate(0b00,0b00) = 1;
    gate(0b01,0b10) = 1;
    gate(0b10,0b01) = 1;
    gate(0b11,0b11) = 1;

    // Apply SWAP to qubits ctrl = 1, targ = 0
    unsigned ctrl = 1;
    unsigned targ = 0;
    const auto m0 = makeMatrix(gate, 3, ctrl, targ);
    std::cout << m0 << std::endl;
    // Check all the basis states
    const std::complex<double> one{1,0};
    EXPECT_EQ(m0(0b000,0b000), one); // |000) -> |000) 
    EXPECT_EQ(m0(0b001,0b010), one); // |001) -> |010) 
    EXPECT_EQ(m0(0b010,0b001), one); // |010) -> |001) 
    EXPECT_EQ(m0(0b011,0b011), one); // |011) -> |011) 
    EXPECT_EQ(m0(0b100,0b100), one); // |100) -> |100) 
    EXPECT_EQ(m0(0b101,0b110), one); // |101) -> |110) 
    EXPECT_EQ(m0(0b110,0b101), one); // |110) -> |101) 
    EXPECT_EQ(m0(0b111,0b111), one); // |111) -> |111) 

    // Apply SWAP to qubits ctrl = 0, targ = 2
    ctrl = 0;
    targ = 2;
    const auto m1 = makeMatrix(gate, 3, ctrl, targ);
    std::cout << m1 << std::endl;
    // Check all the basis states
    EXPECT_EQ(m1(0b000,0b000), one);
    EXPECT_EQ(m1(0b001,0b100), one);
    EXPECT_EQ(m1(0b010,0b010), one);
    EXPECT_EQ(m1(0b011,0b110), one);
    EXPECT_EQ(m1(0b100,0b001), one);
    EXPECT_EQ(m1(0b101,0b101), one);
    EXPECT_EQ(m1(0b110,0b011), one);
    EXPECT_EQ(m1(0b111,0b111), one);

}

TEST(GateTests,MakeMatrixTestsCHadamard)
{
    // Controlled Hadamard gate, qubit 1 is ctrl, qubit 0 is targ
    arma::Mat<std::complex<double>> gate(4, 4, arma::fill::zeros);
    double one_sqrt2{1/std::sqrt(2)};
    gate(0b00,0b00) = 1;
    gate(0b01,0b01) = 1;
    gate(0b10,0b10) = one_sqrt2;
    gate(0b10,0b11) = one_sqrt2;
    gate(0b11,0b10) = one_sqrt2;
    gate(0b11,0b11) = -one_sqrt2;

    // Apply controlled Hadamard to qubits ctrl = 1, targ = 0
    unsigned ctrl = 1;
    unsigned targ = 0;
    const auto m0 = makeMatrix(gate, 3, ctrl, targ);
    std::cout << m0 << std::endl;
    // Check all the basis states
    const std::complex<double> one{1,0};
    const std::complex<double> val{one_sqrt2,0};
    EXPECT_EQ(m0(0b000,0b000), one); // |000) -> |000) 
    EXPECT_EQ(m0(0b001,0b001), one); // |001) -> |010) 
    EXPECT_EQ(m0(0b010,0b010), val); // |010) -> [ |010) + |011) ] / sqrt2 
    EXPECT_EQ(m0(0b010,0b011), val);
    EXPECT_EQ(m0(0b011,0b010), val); // |011) -> [ |010) - |011) ] / sqrt2 
    EXPECT_EQ(m0(0b011,0b011), -val);
    EXPECT_EQ(m0(0b100,0b100), one); // |100) -> |100) 
    EXPECT_EQ(m0(0b101,0b101), one); // |101) -> |110) 
    EXPECT_EQ(m0(0b110,0b110), val); // |110) -> [ |110) + |111) ] / sqrt2 
    EXPECT_EQ(m0(0b110,0b111), val);
    EXPECT_EQ(m0(0b111,0b110), val); // |111) -> [ |110) - |111) ] / sqrt2 
    EXPECT_EQ(m0(0b111,0b111), -val);

    // Apply controlled Hadamard to qubits ctrl = 1, targ = 0
    ctrl = 0;
    targ = 2;
    const auto m1 = makeMatrix(gate, 3, ctrl, targ);
    std::cout << m1 << std::endl;

    // Check all the basis states
    // Identity
    EXPECT_EQ(m1(0b000,0b000), one);
    EXPECT_EQ(m1(0b010,0b010), one);
    EXPECT_EQ(m1(0b100,0b100), one);
    EXPECT_EQ(m1(0b110,0b110), one);

    // To the plus state
    EXPECT_EQ(m1(0b001,0b001), val);
    EXPECT_EQ(m1(0b001,0b101), val);
    EXPECT_EQ(m1(0b011,0b011), val);
    EXPECT_EQ(m1(0b011,0b111), val);

    // To the minus state
    EXPECT_EQ(m1(0b101,0b001), val);
    EXPECT_EQ(m1(0b101,0b101), -val);
    EXPECT_EQ(m1(0b111,0b011), val);
    EXPECT_EQ(m1(0b111,0b111), -val);

}

/// Test the makeMatrix function against tensor product matrix
TEST(GateTests,MakeMatrixTestsTensor)
{
    const unsigned num_qubits{ 8 };
    unsigned ctrl = 6;
    unsigned targ = 5;
    
    // Controlled Hadamard gate, qubit 1 is ctrl, qubit 0 is targ
    arma::Mat<std::complex<double>> mat(4, 4, arma::fill::zeros);
    double one_sqrt2{1/std::sqrt(2)};
    mat(0b00,0b00) = 1;
    mat(0b01,0b01) = 1;
    mat(0b10,0b10) = one_sqrt2;
    mat(0b10,0b11) = one_sqrt2;
    mat(0b11,0b10) = one_sqrt2;
    mat(0b11,0b11) = -one_sqrt2;

    // Make identity gate
    // Create gate in armadillo
    // Sizes of idenity matrix padding
    std::size_t pre = 1 << targ;
    std::size_t post = 1 << (num_qubits - ctrl - 1);

    // Tensor to make the gate. Have fun with these lines...
    arma::SpMat<std::complex<double>> gate =
	arma::speye<arma::SpMat<std::complex<double>>>(pre, pre);
    gate = arma::kron(arma::conv_to<arma::SpMat<std::complex<double>>>::from(mat),
		      gate);
    gate = arma::kron(arma::speye<arma::SpMat<std::complex<double>>>(post, post),
		      gate);


    // Compute the same matrix using the makeMatrix function
    const auto gate2{ makeMatrix(mat, num_qubits, ctrl, targ) };
    EXPECT_TRUE(arma::approx_equal(gate2, gate, "absdiff", 1e-8));    
}



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

TYPED_TEST(NPOneQubitGates, OneQubitNoArg)
{   
    const unsigned num_qubits = 8;
    const unsigned targ = 4;
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;

    // Create list of gates mapped to matrices
    std::vector<std::pair<std::function<void(Sim &, unsigned)>, 
			  arma::SpMat<std::complex<Fp>>>> gates;

    // PauliZ
    arma::SpMat<std::complex<Fp>> pauliZ(2, 2);
    pauliZ(0, 0) = 1;
    pauliZ(1, 1) = -1;
    auto fn_pauliZ = [](Sim & sim, unsigned targ) {
			 sim.pauliZ(targ);
		     };
    gates.push_back({fn_pauliZ, pauliZ});
    
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
