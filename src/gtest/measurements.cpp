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

#define UNCERTAINTY 0.1

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
	
	    // Measure, get outcome, and collapse the state
	    unsigned outcome = q.measure(n);

	    // Update the statistics
	    X(outcome);
	    
	    // Check that the state collapsed to the correct thing
	    arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	    arma::Col<std::complex<Fp>> state_2{ applyProjector(P[outcome], v) };
	    EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));
	}

	// Compute the probability of getting 1
	Fp p1_arma = probability(P[1],v);
	Fp p1 = X.mean();
	
	// Compare with the true mean. This test will succeed if the
	// estimated probability is within 5% of the true value.
	///\todo We need to figure out a legitimate way to test
	/// whether the probability is correct.
	EXPECT_NEAR(p1, p1_arma, UNCERTAINTY);		
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
	
	    // Measure, get outcome, and collapse the state
	    unsigned outcome = q.measure(n);

	    // Update the statistics
	    X(outcome);
	    
	    // Check that the state collapsed to the correct thing
	    arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	    arma::Col<std::complex<Fp>> state_2{ applyProjector(P[outcome], v) };
	    EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));
	}

	// Compute the probability of getting 1
	Fp p1_arma = probability(P[1],v);
	Fp p1 = X.mean();
	
	// Compare with the true mean. This test will succeed if the
	// estimated probability is within 5% of the true value.
	///\todo We need to figure out a legitimate way to test
	/// whether the probability is correct.
	EXPECT_NEAR(p1, p1_arma, UNCERTAINTY);	
	
	
    }
}

TYPED_TEST(Measurements, SampleTest)
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

    for (unsigned n = 0; n < num_qubits; n++) {

	// Make the projector onto the |0) state
	arma::SpMat<std::complex<Fp>> P0 { projector<Fp>(num_qubits, n, 0) };

	// Sample the measure function many times
	std::size_t samples{ 1000 };

	// Sample outcomes from the state vector
	std::vector<std::size_t> outcomes = q.sample(n, samples);

	// Compute the probability of getting 0
	Fp p0 = probability(P0,v);
	Fp p0_arma = static_cast<double>(outcomes[0])/(outcomes[0] + outcomes[1]);
	
	// Compare with the true mean. This test will succeed if the
	// estimated probability is within 5% of the true value.
	///\todo We need to figure out a legitimate way to test
	/// whether the probability is correct.
	EXPECT_NEAR(p0, p0_arma, UNCERTAINTY);	
    }
}

TYPED_TEST(NPMeasurements, SampleTest)
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

    for (unsigned n = 0; n < num_qubits; n++) {

	// Make the projector onto the |0) state
	arma::SpMat<std::complex<Fp>> P0 { projector<Fp>(num_qubits, n, 0) };

	// Sample the measure function many times
	std::size_t samples{ 1000 };

	// Sample outcomes from the state vector
	std::vector<std::size_t> outcomes = q.sample(n, samples);

	// Compute the probability of getting 0
	Fp p0 = probability(P0,v);
	Fp p0_arma = static_cast<double>(outcomes[0])/(outcomes[0] + outcomes[1]);
	
	// Compare with the true mean. This test will succeed if the
	// estimated probability is within 5% of the true value.
	///\todo We need to figure out a legitimate way to test
	/// whether the probability is correct.
	EXPECT_NEAR(p0, p0_arma, UNCERTAINTY);	
    }
}

TYPED_TEST(Measurements, MeasureAllTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 4 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomState<Fp>(num_qubits);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };
	
    // Sample the measure function many times
    std::size_t samples{ 1000 };

    // Record how many times each outcome was measured
    std::map<std::size_t, unsigned> results;
    
    for (std::size_t s = 0; s < samples; s++) {

	// Reset the state to the copy
	q = q_copy;
	
	// Measure, get outcome bitstring, and collapse to basis state
	std::size_t outcome = q.measureAll();
	results[outcome]++;

	// Make the resulting basis state
	const std::size_t dim{ 1ULL << num_qubits  };
	arma::Col<std::complex<Fp>> state_2{ dim, arma::fill::zeros };
	state_2(outcome) = 1;
	
	// Check that the state collapsed to the correct thing
	arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));
    }

    // Check that the probability of each outcome is correct
    for (const auto & [outcome, s] : results) {

	// Only check stats if there are enough samples
	if (s > 50) {

	    // Compute the probability estimate
	    double p_est = static_cast<double>(s)/samples;

	    // Compute the true probability
	    double p_true = std::abs(v[outcome]) * std::abs(v[outcome]);

	    // Check that the estimate is near the true value (with 5%)
	    ///\todo Find a legitimate statistical test.
	    EXPECT_NEAR(p_est, p_true, UNCERTAINTY);

	}

    }
}


TYPED_TEST(NPMeasurements, MeasureAllTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 4 };
    const unsigned num_ones{ 2 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomNPState<Fp>(num_qubits, num_ones);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };
	
    // Sample the measure function many times
    std::size_t samples{ 1000 };

    // Record how many times each outcome was measured
    std::map<std::size_t, unsigned> results;
    
    for (std::size_t s = 0; s < samples; s++) {

	// Reset the state to the copy
	q = q_copy;
	
	// Measure, get outcome bitstring, and collapse to basis state
	std::size_t outcome = q.measureAll();
	results[outcome]++;

	// Make the resulting basis state
	const std::size_t dim{ 1ULL << num_qubits  };
	arma::Col<std::complex<Fp>> state_2{ dim, arma::fill::zeros };
	state_2(outcome) = 1;
	
	// Check that the state collapsed to the correct thing
	arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));
    }

    // Check that the probability of each outcome is correct
    for (const auto & [outcome, s] : results) {

	// Only check stats if there are enough samples
	if (s > 50) {

	    // Compute the probability estimate
	    double p_est = static_cast<double>(s)/samples;

	    // Compute the true probability
	    double p_true = std::abs(v[outcome]) * std::abs(v[outcome]);

	    // Check that the estimate is near the true value (with 5%)
	    ///\todo Find a legitimate statistical test.
	    EXPECT_NEAR(p_est, p_true, UNCERTAINTY);

	}

    }
}

TYPED_TEST(Measurements, SampleAllTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 4 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomState<Fp>(num_qubits);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };
	
    // Sample the measure function many times
    std::size_t samples{ 1000 };

    // Record how many times each outcome was measured
    std::map<std::size_t, std::size_t> results = q.sampleAll(samples);
    
    // Check that the probability of each outcome is correct
    for (const auto & [outcome, s] : results) {

	// Only check stats if there are enough samples
	if (s > 50) {

	    // Compute the probability estimate
	    double p_est = static_cast<double>(s)/samples;

	    // Compute the true probability
	    double p_true = std::abs(v[outcome]) * std::abs(v[outcome]);

	    // Check that the estimate is near the true value (with 5%)
	    ///\todo Find a legitimate statistical test.
	    EXPECT_NEAR(p_est, p_true, UNCERTAINTY);

	}

    }
}

TYPED_TEST(NPMeasurements, SampleAllTest)
{
    using Sim = TypeParam::Sim;
    using Fp = TypeParam::Sim::Fp_type;
    const unsigned num_qubits{ 4 };
    const unsigned num_ones{ 2 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomNPState<Fp>(num_qubits, num_ones);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };
	
    // Sample the measure function many times
    std::size_t samples{ 1000 };

    // Record how many times each outcome was measured
    std::map<std::size_t, std::size_t> results = q.sampleAll(samples);
    
    // Check that the probability of each outcome is correct
    for (const auto & [outcome, s] : results) {

	// Only check stats if there are enough samples
	if (s > 50) {

	    // Compute the probability estimate
	    double p_est = static_cast<double>(s)/samples;

	    // Compute the true probability
	    double p_true = std::abs(v[outcome]) * std::abs(v[outcome]);

	    // Check that the estimate is near the true value (with 5%)
	    ///\todo Find a legitimate statistical test.
	    EXPECT_NEAR(p_est, p_true, UNCERTAINTY);

	}

    }
}

TEST(SampleAll2Test, DefaultDouble)
{
    using Sim = qsl::Qubits<qsl::Type::Default, double>;
    using Fp = Sim::Fp_type;
    const unsigned num_qubits{ 4 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomState<Fp>(num_qubits);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };
	
    // Sample the measure function many times
    std::size_t samples{ 1000 };

    // Record how many times each outcome was measured
    std::map<std::size_t, std::size_t> results = q.sampleAll2(samples);
    
    // Check that the probability of each outcome is correct
    for (const auto & [outcome, s] : results) {

	// Only check stats if there are enough samples
	if (s > 50) {

	    // Compute the probability estimate
	    double p_est = static_cast<double>(s)/samples;

	    // Compute the true probability
	    double p_true = std::abs(v[outcome]) * std::abs(v[outcome]);

	    // Check that the estimate is near the true value (with 5%)
	    ///\todo Find a legitimate statistical test.
	    EXPECT_NEAR(p_est, p_true, UNCERTAINTY);

	}

   } 
}

TEST(SampleAll2Test, ResizeDouble)
{
    using Sim = qsl::Qubits<qsl::Type::Resize, double>;
    using Fp = Sim::Fp_type;
    const unsigned num_qubits{ 4 };

    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomState<Fp>(num_qubits);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };
	
    // Sample the measure function many times
    std::size_t samples{ 1000 };

    // Record how many times each outcome was measured
    std::map<std::size_t, std::size_t> results = q.sampleAll2(samples);
    
    // Check that the probability of each outcome is correct
    for (const auto & [outcome, s] : results) {

	// Only check stats if there are enough samples
	if (s > 50) {

	    // Compute the probability estimate
	    double p_est = static_cast<double>(s)/samples;

	    // Compute the true probability
	    double p_true = std::abs(v[outcome]) * std::abs(v[outcome]);

	    // Check that the estimate is near the true value (with 5%)
	    ///\todo Find a legitimate statistical test.
	    EXPECT_NEAR(p_est, p_true, UNCERTAINTY);

	}

    }
}


TEST(SampleAll2Test, NPDouble)
{
    using Sim = qsl::Qubits<qsl::Type::NP, double>;
    using Fp = Sim::Fp_type;
    const unsigned num_qubits{ 4 };
    const unsigned num_ones{ 2 };
    // Make a random state
    Sim q{num_qubits};
    const std::vector<qsl::complex<Fp>> state
	= qsl::makeRandomNPState<Fp>(num_qubits, num_ones);
    q.setState(state);

    // Make a copy to re-initialise the state
    const Sim q_copy{q}; 
    
    // Set an armadillo vector to the same state
    arma::Col<std::complex<Fp>> v{ toArmaState(q) };
	
    // Sample the measure function many times
    std::size_t samples{ 1000 };

    // Record how many times each outcome was measured
    std::map<std::size_t, std::size_t> results = q.sampleAll2(samples);
    
    // Check that the probability of each outcome is correct
    for (const auto & [outcome, s] : results) {

	// Only check stats if there are enough samples
	if (s > 50) {

	    // Compute the probability estimate
	    double p_est = static_cast<double>(s)/samples;

	    // Compute the true probability
	    double p_true = std::abs(v[outcome]) * std::abs(v[outcome]);

	    // Check that the estimate is near the true value (with 5%)
	    ///\todo Find a legitimate statistical test.
	    EXPECT_NEAR(p_est, p_true, UNCERTAINTY);

	}

    }
}

