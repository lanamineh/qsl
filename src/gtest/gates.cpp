#include <gtest/gtest.h>
#include <qsl/qubits.hpp>
#include <qsl/concepts.hpp>
#include <armadillo>
#include <complex>
#include <list>
#include "test-utils.hpp"

template<typename T>
struct SimWrapper
{
    using Sim = T;
};

/**
 * \brief Typed test suite for one-qubit gates
 */
template <typename T>
class Gates : public testing::Test {
public:
    using List = std::list<T>;
    static T shared_;
    T value_;
};

template <typename T>
class NPGates : public testing::Test {
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

TYPED_TEST_SUITE(Gates, SimTypes);
TYPED_TEST_SUITE(NPGates, NPSimTypes);

TEST(GateTests, MakeMatrixTestsPauliX)
{
    // Pauli gate, targ = 4
    const unsigned num_qubits{ 7 };
    const unsigned targ{ 4 };
    arma::Mat<std::complex<double>> mat(2, 2, arma::fill::zeros);
    mat(0b0,0b1) = 1;
    mat(0b1,0b0) = 1;
    
    // Make identity gate
    // Create gate in armadillo
    // Sizes of idenity matrix padding
    std::size_t pre = 1 << targ;
    std::size_t post = 1 << (num_qubits - targ - 1);

    // Tensor to make the gate. Have fun with these lines...
    arma::SpMat<std::complex<double>> gate =
	arma::speye<arma::SpMat<std::complex<double>>>(pre, pre);
    gate = arma::kron(arma::conv_to<arma::SpMat<std::complex<double>>>::from(mat),
		      gate);
    gate = arma::kron(arma::speye<arma::SpMat<std::complex<double>>>(post, post),
		      gate);

    // Compute the same matrix using the makeMatrix function
    const auto gate2{ makeMatrix(mat, num_qubits, {targ}) };
    EXPECT_TRUE(arma::approx_equal(gate2, gate, "absdiff", 1e-8));    



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
    const auto m0 = makeMatrix(gate, 3, {targ, ctrl});
    //std::cout << m0 << std::endl;
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
    const auto m1 = makeMatrix(gate, 3, {targ, ctrl});
    //std::cout << m1 << std::endl;
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
    const auto m2 = makeMatrix(gate, 2, {targ, ctrl});
    //std::cout << m2 << std::endl;
    // Check all the basis states
    EXPECT_EQ(m2(0b00,0b00), one); 
    EXPECT_EQ(m2(0b01,0b11), one);
    EXPECT_EQ(m2(0b10,0b10), one);
    EXPECT_EQ(m2(0b11,0b01), one);

    // Check the degnerate case the other way round
    ctrl = 1;
    targ = 0;
    const auto m3 = makeMatrix(gate, 2, {targ, ctrl});
    //std::cout << m3 << std::endl;
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
    const auto m0 = makeMatrix(gate, 3, {targ, ctrl});
    //std::cout << m0 << std::endl;
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
    const auto m1 = makeMatrix(gate, 3, {targ, ctrl});
    //std::cout << m1 << std::endl;
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
    const auto m0 = makeMatrix(gate, 3, {targ, ctrl});
    //std::cout << m0 << std::endl;
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
    const auto m1 = makeMatrix(gate, 3, {targ, ctrl});
    //std::cout << m1 << std::endl;

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
    const auto gate2{ makeMatrix(mat, num_qubits, {targ, ctrl}) };
    EXPECT_TRUE(arma::approx_equal(gate2, gate, "absdiff", 1e-8));    
}



TYPED_TEST(Gates, OneQubitGate)
{   
    const unsigned num_qubits = 8;
    const unsigned targ = 4;
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    // Generate random angle to test parameterised gates
    Fp angle = arma::datum::pi * arma::randn<Fp>();
    
    // Create list of gates mapped to matrices
    std::vector<std::pair<std::function<void(Sim &, unsigned)>, 
			  arma::Mat<std::complex<Fp>>>> gates;

    // PauliX 
    arma::Mat<std::complex<Fp>> pauliX(2, 2, arma::fill::zeros);
    pauliX(0, 1) = 1;
    pauliX(1, 0) = 1;
    auto fn_pauliX = [](Sim & sim, unsigned targ) {
			 sim.pauliX(targ);
		     };
    gates.push_back({fn_pauliX, pauliX});

    // PauliY 
    arma::Mat<std::complex<Fp>> pauliY(2, 2, arma::fill::zeros);
    pauliY(0, 1) = std::complex<Fp>{0, -1};
    pauliY(1, 0) = std::complex<Fp>{0, 1};
    auto fn_pauliY = [](Sim & sim, unsigned targ) {
    			 sim.pauliY(targ);
    		     };
    gates.push_back({fn_pauliY, pauliY});

    // PauliZ
    arma::Mat<std::complex<Fp>> pauliZ(2, 2, arma::fill::zeros);
    pauliZ(0, 0) = 1;
    pauliZ(1, 1) = -1;
    auto fn_pauliZ = [](Sim & sim, unsigned targ) {
			 sim.pauliZ(targ);
		     };
    gates.push_back({fn_pauliZ, pauliZ});
    
    // Hadamard
    arma::Mat<std::complex<Fp>> hadamard(2, 2, arma::fill::zeros);
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
    arma::Mat<std::complex<Fp>> phase(2, 2, arma::fill::zeros);
    phase(0, 0) = 1;
    phase(1, 1) = std::complex<Fp>{std::cos(angle), std::sin(angle)};
    auto fn_phase = [=](Sim & sim, unsigned targ) {
			sim.phase(targ, angle);
		       };
    gates.push_back({fn_phase, phase});
    
    // rotateX
    arma::Mat<std::complex<Fp>> rotateX(2, 2, arma::fill::zeros);
    rotateX(0, 0) = std::cos(angle/2);
    rotateX(0, 1) = std::complex<Fp>{0,-std::sin(angle/2)};
    rotateX(1, 0) = std::complex<Fp>{0,-std::sin(angle/2)};
    rotateX(1, 1) = std::cos(angle/2);
    auto fn_rotateX = [=](Sim & sim, unsigned targ) {
			sim.rotateX(targ, angle);
		       };
    gates.push_back({fn_rotateX, rotateX});

    // rotateY
    arma::Mat<std::complex<Fp>> rotateY(2, 2, arma::fill::zeros);
    rotateY(0, 0) = std::cos(angle/2);
    rotateY(0, 1) = -std::sin(angle/2);
    rotateY(1, 0) = std::sin(angle/2);
    rotateY(1, 1) = std::cos(angle/2);
    auto fn_rotateY = [=](Sim & sim, unsigned targ) {
			sim.rotateY(targ, angle);
		       };
    gates.push_back({fn_rotateY, rotateY});
    
    // rotateZ    
    arma::Mat<std::complex<Fp>> rotateZ(2, 2, arma::fill::zeros);
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

	// Apply gate in armadillo
	v = makeMatrix(mat, num_qubits, {targ}) * v;

	// Read qubit state into armadillo
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}
    
    	EXPECT_TRUE(arma::approx_equal(v, qubit_v, "both", 1e-6, 1e-8));
    }
}

TYPED_TEST(Gates, TwoQubitGate)
{   
    const unsigned num_qubits = 8;
    const unsigned ctrl = 2;
    const unsigned targ = 4;
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    // Generate random angle to test parameterised gates
    Fp angle = arma::datum::pi * arma::randn<Fp>();
    
    // Create list of gates mapped to matrices
    std::vector<std::pair<std::function<void(Sim &, unsigned, unsigned)>, 
			  arma::Mat<std::complex<Fp>>>> gates;

    // controlNot
    arma::Mat<std::complex<Fp>> cnot(4, 4, arma::fill::zeros);
    cnot(0b00,0b00) = 1;
    cnot(0b01,0b01) = 1;
    cnot(0b10,0b11) = 1;
    cnot(0b11,0b10) = 1;
    auto fn_cnot = [](Sim & sim, unsigned ctrl, unsigned targ) {
			 sim.controlNot(ctrl, targ);
		     };
    gates.push_back({fn_cnot, cnot});

    // controlY
    arma::Mat<std::complex<Fp>> cy(4, 4, arma::fill::zeros);
    cy(0b00,0b00) = 1;
    cy(0b01,0b01) = 1;
    cy(0b10,0b11) = std::complex<Fp>{0, -1};
    cy(0b11,0b10) = std::complex<Fp>{0, 1};;
    auto fn_cy = [=](Sim & sim, unsigned ctrl, unsigned targ) {
		     sim.controlY(ctrl, targ);
		   };
    gates.push_back({fn_cy, cy});

    // controlZ
    arma::Mat<std::complex<Fp>> cz(4, 4, arma::fill::zeros);
    cz(0b00,0b00) = 1;
    cz(0b01,0b01) = 1;
    cz(0b10,0b10) = 1;
    cz(0b11,0b11) = -1;
    auto fn_cz = [=](Sim & sim, unsigned ctrl, unsigned targ) {
		     sim.controlZ(ctrl, targ);
		   };
    gates.push_back({fn_cz, cz});

    // controlRotateX
    arma::Mat<std::complex<Fp>> crx(4, 4, arma::fill::zeros);
    crx(0b00,0b00) = 1;
    crx(0b01,0b01) = 1;
    crx(0b10,0b10) = std::complex<Fp>{std::cos(angle/2), 0};
    crx(0b10,0b11) = std::complex<Fp>{0, -std::sin(angle/2)};
    crx(0b11,0b10) = std::complex<Fp>{0, -std::sin(angle/2)};
    crx(0b11,0b11) = std::complex<Fp>{std::cos(angle/2), 0};

    auto fn_crx = [=](Sim & sim, unsigned ctrl, unsigned targ) {
	sim.controlRotateX(ctrl, targ, angle);
		   };
    gates.push_back({fn_crx, crx});

    // controlRotateY
    arma::Mat<std::complex<Fp>> cry(4, 4, arma::fill::zeros);
    cry(0b00,0b00) = 1;
    cry(0b01,0b01) = 1;
    cry(0b10,0b10) = std::cos(angle/2);
    cry(0b10,0b11) = -std::sin(angle/2);
    cry(0b11,0b10) = std::sin(angle/2);
    cry(0b11,0b11) = std::cos(angle/2);

    auto fn_cry = [=](Sim & sim, unsigned ctrl, unsigned targ) {
	sim.controlRotateY(ctrl, targ, angle);
		   };
    gates.push_back({fn_cry, cry});

    // controlRotateZ
    arma::Mat<std::complex<Fp>> crz(4, 4, arma::fill::zeros);
    crz(0b00,0b00) = 1;
    crz(0b01,0b01) = 1;
    crz(0b10,0b10) = std::complex<Fp>{std::cos(angle/2), -std::sin(angle/2)};
    crz(0b11,0b11) = std::complex<Fp>{std::cos(angle/2), std::sin(angle/2)};
    auto fn_crz = [=](Sim & sim, unsigned ctrl, unsigned targ) {
	sim.controlRotateZ(ctrl, targ, angle);
		   };
    gates.push_back({fn_crz, crz});
    
    // controlPhase
    arma::Mat<std::complex<Fp>> cphase(4, 4, arma::fill::zeros);
    cphase(0b00,0b00) = 1;
    cphase(0b01,0b01) = 1;
    cphase(0b10,0b10) = 1;
    cphase(0b11,0b11) = std::complex<Fp>{std::cos(angle), std::sin(angle)};
    auto fn_cphase = [=](Sim & sim, unsigned ctrl, unsigned targ) {
		       sim.controlPhase(ctrl, targ, angle);
		   };
    gates.push_back({fn_cphase, cphase});

    // controlHadamard
    Fp sqrt2 = 1/std::sqrt(2);
    arma::Mat<std::complex<Fp>> ch(4, 4, arma::fill::zeros);
    ch(0b00,0b00) = 1;
    ch(0b01,0b01) = 1;
    ch(0b10,0b10) = sqrt2;
    ch(0b10,0b11) = sqrt2;
    ch(0b11,0b10) = sqrt2;
    ch(0b11,0b11) = -sqrt2;

    auto fn_ch = [=](Sim & sim, unsigned ctrl, unsigned targ) {
	sim.controlHadamard(ctrl, targ);
		   };
    gates.push_back({fn_ch, ch});
    
    // swap
    arma::Mat<std::complex<Fp>> swap(4, 4, arma::fill::zeros);
    swap(0b00,0b00) = 1;
    swap(0b01,0b10) = 1;
    swap(0b10,0b01) = 1;
    swap(0b11,0b11) = 1;
    auto fn_swap = [=](Sim & sim, unsigned ctrl, unsigned targ) {
		       sim.swap(ctrl, targ);
		   };
    gates.push_back({fn_swap, swap});

    // fswap
    arma::Mat<std::complex<Fp>> fswap(4, 4, arma::fill::zeros);
    fswap(0b00,0b00) = 1;
    fswap(0b01,0b10) = 1;
    fswap(0b10,0b01) = 1;
    fswap(0b11,0b11) = -1;
    auto fn_fswap = [=](Sim & sim, unsigned ctrl, unsigned targ) {
	sim.fswap(ctrl, targ);
		   };
    gates.push_back({fn_fswap, fswap});

    
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
	std::invoke(fn, q, ctrl, targ);

	// Apply gate in armadillo
	v = makeMatrix(mat, num_qubits, {targ, ctrl}) * v;

	// Read qubit state into armadillo
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}
    
    	EXPECT_TRUE(arma::approx_equal(v, qubit_v, "both", 1e-6, 1e-8));
    }
}

TYPED_TEST(NPGates, TwoQubitGate)
{   
    const unsigned num_qubits = 8;
    const unsigned num_ones = 4;
    const unsigned ctrl = 2;
    const unsigned targ = 4;
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    // Generate random angle to test parameterised gates
    Fp angle = arma::datum::pi * arma::randn<Fp>();
    
    // Create list of gates mapped to matrices
    std::vector<std::pair<std::function<void(Sim &, unsigned, unsigned)>, 
			  arma::Mat<std::complex<Fp>>>> gates;

    // controlZ
    arma::Mat<std::complex<Fp>> cz(4, 4, arma::fill::zeros);
    cz(0b00,0b00) = 1;
    cz(0b01,0b01) = 1;
    cz(0b10,0b10) = 1;
    cz(0b11,0b11) = -1;
    auto fn_cz = [=](Sim & sim, unsigned ctrl, unsigned targ) {
		     sim.controlZ(ctrl, targ);
		   };
    gates.push_back({fn_cz, cz});

    // controlRotateZ
    arma::Mat<std::complex<Fp>> crz(4, 4, arma::fill::zeros);
    crz(0b00,0b00) = 1;
    crz(0b01,0b01) = 1;
    crz(0b10,0b10) = std::complex<Fp>{std::cos(angle/2), -std::sin(angle/2)};
    crz(0b11,0b11) = std::complex<Fp>{std::cos(angle/2), std::sin(angle/2)};
    auto fn_crz = [=](Sim & sim, unsigned ctrl, unsigned targ) {
	sim.controlRotateZ(ctrl, targ, angle);
		   };
    gates.push_back({fn_crz, crz});
    
    // controlPhase
    arma::Mat<std::complex<Fp>> cphase(4, 4, arma::fill::zeros);
    cphase(0b00,0b00) = 1;
    cphase(0b01,0b01) = 1;
    cphase(0b10,0b10) = 1;
    cphase(0b11,0b11) = std::complex<Fp>{std::cos(angle), std::sin(angle)};
    auto fn_cphase = [=](Sim & sim, unsigned ctrl, unsigned targ) {
		       sim.controlPhase(ctrl, targ, angle);
		   };
    gates.push_back({fn_cphase, cphase});

    // swap
    arma::Mat<std::complex<Fp>> swap(4, 4, arma::fill::zeros);
    swap(0b00,0b00) = 1;
    swap(0b01,0b10) = 1;
    swap(0b10,0b01) = 1;
    swap(0b11,0b11) = 1;
    auto fn_swap = [=](Sim & sim, unsigned ctrl, unsigned targ) {
		       sim.swap(ctrl, targ);
		   };
    gates.push_back({fn_swap, swap});

    // fswap
    arma::Mat<std::complex<Fp>> fswap(4, 4, arma::fill::zeros);
    fswap(0b00,0b00) = 1;
    fswap(0b01,0b10) = 1;
    fswap(0b10,0b01) = 1;
    fswap(0b11,0b11) = -1;
    auto fn_fswap = [=](Sim & sim, unsigned ctrl, unsigned targ) {
	sim.fswap(ctrl, targ);
		   };
    gates.push_back({fn_fswap, fswap});

    
    for (const auto & [fn, mat] : gates) {    

	// Make a random state
	Sim q{num_qubits};
	const std::vector<qsl::complex<Fp>> state
	    = qsl::makeRandomNPState<Fp>(num_qubits, num_ones);
	q.setState(state);

	// Set an armadillo vector to the same state
	std::size_t dim = 1 << num_qubits;
	arma::Col<std::complex<Fp>> v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    v(i) = std::complex<Fp>{state[i].real, state[i].imag};
	}

	// Apply gate to qubits
	std::invoke(fn, q, ctrl, targ);

	// Apply gate in armadillo
	v = makeMatrix(mat, num_qubits, {targ, ctrl}) * v;

	// Read qubit state into armadillo
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}
    
    	EXPECT_TRUE(arma::approx_equal(v, qubit_v, "both", 1e-6, 1e-8));
    }
}


TYPED_TEST(NPGates, OneQubitGate)
{   
    const unsigned num_qubits = 8;
    const unsigned num_ones = 4;
    const unsigned targ = 4;
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    // Generate random angle to test parameterised gates
    Fp angle = arma::datum::pi * arma::randn<Fp>();
    
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
	    = qsl::makeRandomNPState<Fp>(num_qubits, num_ones);
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
