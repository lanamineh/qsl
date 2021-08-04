/*
 *  Copyright 2021 Lana Mineh and John Scott
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */
 
/**
 * \file resize-test.cpp
 * \brief Tests for the resize simulator
 */

#include <gtest/gtest.h>

#include <qsl/qubits.hpp>
#include <qsl/utils.hpp>
#include <list>
#include "test-utils.hpp"

template<typename T>
struct SimWrapper
{
    using Sim = T;
};

/**
 * \brief Typed test suite for simulator basic functions
 *
 * This test suite checks that the simulator basic utility functions
 * work correctly.
 */
template <typename T>
class ResizeTests : public testing::Test {
public:
    using List = std::list<T>;
    static T shared_;
    T value_;
};

/// List of simulator types to check
using Sim1 = SimWrapper<qsl::Qubits<qsl::Type::Resize, float>>;
using Sim2 = SimWrapper<qsl::Qubits<qsl::Type::Resize, double>>;
using SimTypes = ::testing::Types<Sim1, Sim2>;

TYPED_TEST_SUITE(ResizeTests, SimTypes);

/// Qubits object constructors
TYPED_TEST(ResizeTests, ConstructorTest)
{
    using Sim = TypeParam::Sim;
    using Fp_type = TypeParam::Sim::Fp_type;
    
    const unsigned num_qubits{ 5 };
    Sim q{ num_qubits };

    // Set q to a random state
    const std::vector<qsl::complex<Fp_type>> state
	= qsl::makeRandomState<Fp_type>(num_qubits);
    q.setState(state);

    // Construct another object using the copy constructor
    Sim q_copy{ q };    

    // Check that the other properties are also the same
    EXPECT_EQ(q.getNumQubits(), q_copy.getNumQubits());
    
    // Check that both objects have the same internal state
    Fp_type distance = qsl::distance(q.getState(), q_copy.getState());
    EXPECT_FLOAT_EQ(std::abs(distance), 0);

    // Construct another object using initialisation from state vector
    Sim q_from_state{ state };    

    // Check that both objects have the same internal state
    distance = qsl::distance(q.getState(), q_from_state.getState());
    EXPECT_FLOAT_EQ(std::abs(distance), 0);

    // Check that assignment works
    Sim q_new{ num_qubits };
    q_new = q; // Perform the test
    distance = qsl::distance(q.getState(), q_new.getState());
    EXPECT_FLOAT_EQ(std::abs(distance), 0);

    // Check that assignment works for different number of qubits
    Sim q_newer{ num_qubits + 1};
    EXPECT_EQ(q_newer.getNumQubits(), num_qubits + 1); // Perform the test
}

/// Test the basic utility functions of the simulators
TYPED_TEST(ResizeTests, UtilityFunctions)
{
    using Sim = TypeParam::Sim;
    using Fp_type = TypeParam::Sim::Fp_type;
    
    const unsigned num_qubits{ 3 };
    Sim q{ num_qubits };
    
    // Set q to a random state
    const std::vector<qsl::complex<Fp_type>> state
	= qsl::makeRandomState<Fp_type>(num_qubits + 1);

    // Assign state vector with a different number of qubits
    q.setState(state);
    EXPECT_EQ(q.getNumQubits(), num_qubits + 1);
    
    // Check that the basis states work
    std::size_t index = 4; // Should choose a random value
    q.setBasisState(index);
    std::vector<qsl::complex<Fp_type>> basis_state = q.getState();

    // Check that the right element is 1
    EXPECT_FLOAT_EQ(basis_state[index].real, 1);
    EXPECT_FLOAT_EQ(std::abs(norm(basis_state)), 1);
}

/// Check the measure out function of the resize qubit
TYPED_TEST(ResizeTests, MeasureOutTest)
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
	    unsigned outcome = q.measureOut(n);

	    // Update the statistics
	    X(outcome);
	    
	    // Check that the state collapsed to the correct thing
	    arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	    arma::Col<std::complex<Fp>> state_2{ applyProjector(P[outcome], v) };

	    // Remove the nth qubit
	    for (int k = state_2.size()-1; k >= 0; k--) {
		if (getBit(k,n) != outcome) {
		    state_2.shed_row(k);
		}
	    }
		
	    EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));
	}

	// Compute the probability of getting 1
	Fp p1_arma = probability(P[1],v);
	Fp p1 = X.mean();
	
	// Compare with the true mean. This test will succeed if the
	// estimated probability is within 5% of the true value.
	///\todo We need to figure out a legitimate way to test
	/// whether the probability is correct.
	EXPECT_NEAR(p1, p1_arma, 0.05);		
    }
}

TYPED_TEST(ResizeTests, PostselectTest)
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
	Fp p0 = q.postselectOut(n, 0);

	// Calculate the projector for outcome 0
	arma::SpMat<std::complex<Fp>> P0 = projector<Fp>(num_qubits, n, 0);

	// Check that the state collapses to the correct thing
	arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };
	arma::Col<std::complex<Fp>> state_2{ applyProjector(P0, v) };

	// Remove the nth qubit
	for (int k = state_2.size()-1; k >= 0; k--) {
	    if (getBit(k,n) != 0) {
		state_2.shed_row(k);
	    }
	}
	
	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p0_arma = probability(P0,v);	
	EXPECT_NEAR(p0, p0_arma, 1e-6);

	// Reset the state to the copy
	q = q_copy;

	// Postselect on the 1 outcome 
	Fp p1 = q.postselectOut(n, 1);

	// Calculate the projector for outcome 1
	arma::SpMat<std::complex<Fp>> P1 = projector<Fp>(num_qubits, n, 1);

	// Check that the state collapses to the correct thing
	state_1 = toArmaState(q);
	state_2 = applyProjector(P1, v);

	// Remove the nth qubit
	for (int k = state_2.size()-1; k >= 0; k--) {
	    if (getBit(k,n) != 1) {
		state_2.shed_row(k);
	    }
	}

	EXPECT_TRUE(arma::approx_equal(state_1, state_2, "both", 1e-6, 1e-8));

	// Check that the probability returned by postselect is correct
	Fp p1_arma = probability(P1,v);	
	EXPECT_NEAR(p1, p1_arma, 1e-6);	
    }
}

TYPED_TEST(ResizeTests, AddQubitTest)
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

    // Add a qubit in the highest position
    q.addQubit();

    // Check that the state collapses to the correct thing
    arma::Col<std::complex<Fp>> state_1{ toArmaState(q) };

    // Add a qubit by doubling the size of the state vector
    // and filling the second half with zeros
    v.insert_rows(v.size(), v.size(), true);

    // Check for equality
    EXPECT_TRUE(arma::approx_equal(state_1, v, "both", 1e-6, 1e-8));
    
}
