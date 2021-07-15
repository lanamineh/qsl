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
class Measurements : public testing::Test {
public:
    using List = std::list<T>;
    static T shared_;
    T value_;
};

template <typename T>
class NPMeasurements : public testing::Test {
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

TYPED_TEST_SUITE(Measurements, SimTypes);
TYPED_TEST_SUITE(NPMeasurements, NPSimTypes);

TYPED_TEST(Measurements, ProbTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 8 };

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

    // Find the probability of measurement 0 or 1 on every qubit
    for (unsigned n = 0; n < num_qubits; n++) {

	// Find the probability of getting either zero or one
	// from the simulator
	Fp p0 = q.prob(n, 0);
	Fp p1 = q.prob(n, 1);

	// Calculate the projector for the outcome |1), which is
	// the observable that has eigenvalue 0 for |0) and eigenvalue
	// 1 for |1). The expectation value of this observable in a
	// state is the probability of getting 1 on measurement.
	arma::Mat<std::complex<Fp>> projector(2,2,arma::fill::zeros);
	projector(1,1) = 1;
	arma::SpMat<std::complex<Fp>> M = makeMatrix(projector, num_qubits, {n});
	
	// Calculate the probability of zero or one from the
	// armadillo vector
	Fp p1_arma = arma::cdot(v, M*v).real();
	Fp p0_arma = 1 - p1_arma;

	EXPECT_NEAR(p0, p0_arma, 1e-6);
	EXPECT_NEAR(p1, p1_arma, 1e-6);

	// Read qubit state into armadillo to check the state hasn't changed
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}

	EXPECT_TRUE(arma::approx_equal(v, qubit_v, "both", 1e-6, 1e-8));
    }

}

TYPED_TEST(NPMeasurements, ProbTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 8 };
    const unsigned num_ones{ 5 };

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

    // Find the probability of measurement 0 or 1 on every qubit
    for (unsigned n = 0; n < num_qubits; n++) {

	// Find the probability of getting either zero or one
	// from the simulator
	Fp p0 = q.prob(n, 0);
	Fp p1 = q.prob(n, 1);

	// Calculate the projector for the outcome |1), which is
	// the observable that has eigenvalue 0 for |0) and eigenvalue
	// 1 for |1). The expectation value of this observable in a
	// state is the probability of getting 1 on measurement.
	arma::Mat<std::complex<Fp>> projector(2,2,arma::fill::zeros);
	projector(1,1) = 1;
	arma::SpMat<std::complex<Fp>> M = makeMatrix(projector, num_qubits, {n});
	
	// Calculate the probability of zero or one from the
	// armadillo vector
	Fp p1_arma = arma::cdot(v, M*v).real();
	Fp p0_arma = 1 - p1_arma;

	EXPECT_NEAR(p0, p0_arma, 1e-6);
	EXPECT_NEAR(p1, p1_arma, 1e-6);
	
	// Read qubit state into armadillo to check the state hasn't changed
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}

	EXPECT_TRUE(arma::approx_equal(v, qubit_v, "both", 1e-6, 1e-8));
    }

}

TYPED_TEST(Measurements, PostselectTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 8 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomState<Fp>(num_qubits);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    std::size_t dim = 1 << num_qubits;
    arma::Col<std::complex<Fp>> v(dim);
    for (std::size_t i = 0; i < dim; i++) {
	v(i) = std::complex<Fp>{state[i].real, state[i].imag};
    }

    // Find the probability of the 
    for (unsigned n = 0; n < num_qubits; n++) {

	// Reset the state to the copy
	q = q_copy;
	
	// Postselect on the 0 outcome 
	Fp p0 = q.postselect(n, 0);

	// Calculate the projector for outcome 0
	arma::Mat<std::complex<Fp>> proj0(2,2,arma::fill::zeros);
	proj0(0,0) = 1;
	arma::SpMat<std::complex<Fp>> M0 = makeMatrix(proj0, num_qubits, {n});

	// Check that the state collapsed to the right value
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}

	// Check that the state collapses to the correct thing
	arma::Col<std::complex<Fp>> arma_v = (M0 * v)/arma::norm(M0 * v);
	EXPECT_TRUE(arma::approx_equal(arma_v, qubit_v, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p0_arma = arma::cdot(v, M0*v).real();
	EXPECT_NEAR(p0, p0_arma, 1e-6);

	// Reset the state to the copy
	q = q_copy;
	
	// Now postselect on the 1 outcome 
	Fp p1 = q.postselect(n, 1);

	// Calculate the projector for outcome 1
	arma::Mat<std::complex<Fp>> proj1(2,2,arma::fill::zeros);
	proj1(1,1) = 1;
	arma::SpMat<std::complex<Fp>> M1 = makeMatrix(proj1, num_qubits, {n});

	// Check that the state collapsed to the right value
	res = q.getState();
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}

	// Check that the state collapses to the correct thing
        arma_v = (M1 * v)/arma::norm(M1 * v);
	EXPECT_TRUE(arma::approx_equal(arma_v, qubit_v, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p1_arma = arma::cdot(v, M1*v).real();
	EXPECT_NEAR(p1, p1_arma, 1e-6);

	
    }

}

TYPED_TEST(NPMeasurements, PostselectTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 8 };
    const unsigned num_ones{ 5 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomNPState<Fp>(num_qubits, num_ones);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    std::size_t dim = 1 << num_qubits;
    arma::Col<std::complex<Fp>> v(dim);
    for (std::size_t i = 0; i < dim; i++) {
	v(i) = std::complex<Fp>{state[i].real, state[i].imag};
    }

    // Find the probability of the 
    for (unsigned n = 0; n < num_qubits; n++) {

	// Reset the state to the copy
	q = q_copy;
	
	// Postselect on the 0 outcome 
	Fp p0 = q.postselect(n, 0);

	// Calculate the projector for outcome 0
	arma::Mat<std::complex<Fp>> proj0(2,2,arma::fill::zeros);
	proj0(0,0) = 1;
	arma::SpMat<std::complex<Fp>> M0 = makeMatrix(proj0, num_qubits, {n});

	// Check that the state collapsed to the right value
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}

	// Check that the state collapses to the correct thing
	arma::Col<std::complex<Fp>> arma_v = (M0 * v)/arma::norm(M0 * v);
	EXPECT_TRUE(arma::approx_equal(arma_v, qubit_v, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p0_arma = arma::cdot(v, M0*v).real();
	EXPECT_NEAR(p0, p0_arma, 1e-6);

	// Reset the state to the copy
	q = q_copy;
	
	// Now postselect on the 1 outcome 
	Fp p1 = q.postselect(n, 1);

	// Calculate the projector for outcome 1
	arma::Mat<std::complex<Fp>> proj1(2,2,arma::fill::zeros);
	proj1(1,1) = 1;
	arma::SpMat<std::complex<Fp>> M1 = makeMatrix(proj1, num_qubits, {n});

	// Check that the state collapsed to the right value
	res = q.getState();
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}

	// Check that the state collapses to the correct thing
        arma_v = (M1 * v)/arma::norm(M1 * v);
	EXPECT_TRUE(arma::approx_equal(arma_v, qubit_v, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p1_arma = arma::cdot(v, M1*v).real();
	EXPECT_NEAR(p1, p1_arma, 1e-6);

	
    }

}
