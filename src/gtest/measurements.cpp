#include <gtest/gtest.h>
#include <qsl/qubits.hpp>
#include <qsl/concepts.hpp>
#include <armadillo>
#include <complex>
#include <list>
#include "test-utils.hpp"

/**
 * \brief Convert a standard vector state to an armadillo vector
 */
template<std::floating_point Fp>
arma::Col<std::complex<Fp>> toArmaState(const std::vector<qsl::complex<Fp>> & res)
{
    arma::Col<std::complex<Fp>> qubit_v(res.size());
    for (std::size_t i = 0; i < res.size(); i++) {
	qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
    }
    return qubit_v;
}

/**
 * \brief Convert a Simulator state to an armadillo vector
 */
template<qsl::Simulator Sim>
arma::Col<std::complex<typename Sim::Fp_type>> toArmaState(const Sim & sim)
{
    using Fp = Sim::Fp_type;
    // Read qubit state into armadillo to check the state hasn't changed
    std::vector<qsl::complex<Fp>> res = sim.getState();
    return toArmaState<Fp>(res);
}

/**
 * \brief Make the projector onto a particular outcome of the target qubit
 *
 */
template<std::floating_point Fp>
arma::SpMat<std::complex<Fp>>
projector(unsigned num_qubits, unsigned targ, unsigned outcome)
{
    // Calculate the projector for the outcome, which is
    // the observable that has eigenvalue 0 for |~outcome) and eigenvalue
    // 1 for |outcome). The expectation value of this observable in a
    // state is the probability of getting 1 on measurement.
    arma::Mat<std::complex<Fp>> projector(2,2,arma::fill::zeros);
    if (outcome == 0) {
	projector(0,0) = 1;
    } else if (outcome == 1) {
	projector(1,1) = 1;
    } else {
	throw std::out_of_range("outcome must be 0 or 1 in projector() function");
    }
    arma::SpMat<std::complex<Fp>> M = makeMatrix(projector, num_qubits, {targ});
    return M;
}

/**
 * \brief Calculate the probability of measuring a projector outcome
 *
 * Uses the formula prob = (v|P|v), where v is the state and P is 
 * the projector
 */
template<std::floating_point Fp>
Fp probability(const arma::SpMat<std::complex<Fp>> & P,
	       const arma::Col<std::complex<Fp>> & v)
{
    // Calculate the probability of the projector outcome
    Fp prob = arma::cdot(v, P*v).real();
    return prob;
}


/**
 * \brief Collapse a state v using a projector P
 *
 * The output is the state Pv/|Pv| (i.e. the normalised projected state)
 *
 */
template<std::floating_point Fp>
arma::Col<std::complex<Fp>>
applyProjector(const arma::SpMat<std::complex<Fp>> & P,
	       const arma::Col<std::complex<Fp>> & v)
{
    arma::Col<std::complex<Fp>> state = (P * v)/arma::norm(P * v);
    return state;
}

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
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };

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
	arma::SpMat<std::complex<Fp>> P = projector<Fp>(num_qubits, n, 1);
	
	// calculate the probability of zero from the projector P
	Fp p1_arma = probability(P,v);
	Fp p0_arma = 1 - p1_arma;

	EXPECT_NEAR(p0, p0_arma, 1e-6);
	EXPECT_NEAR(p1, p1_arma, 1e-6);

	// Read qubit state into armadillo to check the state hasn't changed
	arma::Col<std::complex<Fp>> qubit_v{ toArmaState(q) };
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
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };

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
	arma::SpMat<std::complex<Fp>> P = projector<Fp>(num_qubits, n, 1);
	
	// calculate the probability of zero from the projector P
	Fp p1_arma = probability(P,v);
	Fp p0_arma = 1 - p1_arma;

	EXPECT_NEAR(p0, p0_arma, 1e-6);
	EXPECT_NEAR(p1, p1_arma, 1e-6);

	// Read qubit state into armadillo to check the state hasn't changed
	arma::Col<std::complex<Fp>> qubit_v{ toArmaState(q) };
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
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };

    for (unsigned n = 0; n < num_qubits; n++) {

	// Reset the state to the copy
	q = q_copy;
	
	// Postselect on the 0 outcome 
	Fp p0 = q.postselect(n, 0);

	// Calculate the projector for outcome 0
	arma::SpMat<std::complex<Fp>> P0 = projector<Fp>(num_qubits, n, 0);

	// Check that the state collapses to the correct thing
	arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	arma::Col<std::complex<Fp>> state_2{ applyProjector(P0, v) };
	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p0_arma = probability(P0,v);	
	EXPECT_NEAR(p0, p0_arma, 1e-6);

	// Reset the state to the copy
	q = q_copy;

	// Postselect on the 1 outcome 
	Fp p1 = q.postselect(n, 1);

	// Calculate the projector for outcome 1
	arma::SpMat<std::complex<Fp>> P1 = projector<Fp>(num_qubits, n, 1);

	// Check that the state collapses to the correct thing
	state_1 = toArmaState(q);
	state_2 = applyProjector(P1, v);
	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p1_arma = probability(P1,v);	
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
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };

    for (unsigned n = 0; n < num_qubits; n++) {

	// Reset the state to the copy
	q = q_copy;
	
	// Postselect on the 0 outcome 
	Fp p0 = q.postselect(n, 0);

	// Calculate the projector for outcome 0
	arma::SpMat<std::complex<Fp>> P0 = projector<Fp>(num_qubits, n, 0);

	// Check that the state collapses to the correct thing
	arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	arma::Col<std::complex<Fp>> state_2{ applyProjector(P0, v) };
	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p0_arma = probability(P0,v);	
	EXPECT_NEAR(p0, p0_arma, 1e-6);

	// Reset the state to the copy
	q = q_copy;

	// Postselect on the 1 outcome 
	Fp p1 = q.postselect(n, 1);

	// Calculate the projector for outcome 1
	arma::SpMat<std::complex<Fp>> P1 = projector<Fp>(num_qubits, n, 1);

	// Check that the state collapses to the correct thing
	state_1 = toArmaState(q);
	state_2 = applyProjector(P1, v);
	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p1_arma = probability(P1,v);	
	EXPECT_NEAR(p1, p1_arma, 1e-6);	
    }
}

TYPED_TEST(Measurements, MeasureTest)
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
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };

    for (unsigned n = 0; n < num_qubits; n++) {

	// Make the projectors for both outcomes
	// The vector is indexed by outcome
	std::vector<arma::SpMat<std::complex<Fp>>> P {
	    projector<Fp>(num_qubits, n, 0),
	    projector<Fp>(num_qubits, n, 1)
	};

	// Store the running average measured outcome
	arma::running_stat<double> X;
	
	// Sample the measure function many times
	std::size_t samples{ 1000 };
	for (std::size_t s = 0; s < samples; s++) {

	    // Reset the state to the copy
	    q = q_copy;
	
	    // Postselect on the 0 outcome 
	    unsigned outcome = q.measure(n);

	    // Update the statistics
	    X(outcome);
	    
	    // Check that the state collapsed to the correct thing
	    arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	    arma::Col<std::complex<Fp>> state_2{ applyProjector(P[outcome], v) };
	    EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));
	}

	// Compute the probability of getting 1
	Fp p1 = probability(P[1],v);
	Fp p1_arma = X.mean();
	
	// Compare with the true mean. This test will succeed if the
	// estimated probability is within 5% of the true value.
	///\todo We need to figure out a legitimate way to test
	/// whether the probability is correct.
	EXPECT_NEAR(p1, p1_arma, 0.05);	
	
	
    }
}

TYPED_TEST(NPMeasurements, MeasureTest)
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
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };

    for (unsigned n = 0; n < num_qubits; n++) {

	// Make the projectors for both outcomes
	// The vector is indexed by outcome
	std::vector<arma::SpMat<std::complex<Fp>>> P {
	    projector<Fp>(num_qubits, n, 0),
	    projector<Fp>(num_qubits, n, 1)
	};

	// Store the running average measured outcome
	arma::running_stat<double> X;
	
	// Sample the measure function many times
	std::size_t samples{ 1000 };
	for (std::size_t s = 0; s < samples; s++) {

	    // Reset the state to the copy
	    q = q_copy;
	
	    // Postselect on the 0 outcome 
	    unsigned outcome = q.measure(n);

	    // Update the statistics
	    X(outcome);
	    
	    // Check that the state collapsed to the correct thing
	    arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	    arma::Col<std::complex<Fp>> state_2{ applyProjector(P[outcome], v) };
	    EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));
	}

	// Compute the probability of getting 1
	Fp p1 = probability(P[1],v);
	Fp p1_arma = X.mean();
	
	// Compare with the true mean. This test will succeed if the
	// estimated probability is within 5% of the true value.
	///\todo We need to figure out a legitimate way to test
	/// whether the probability is correct.
	EXPECT_NEAR(p1, p1_arma, 0.05);	
	
	
    }
}
