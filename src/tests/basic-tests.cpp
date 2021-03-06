/*
 *  Authors: Lana Mineh and John Scott
 *  Copyright 2021 Phasecraft Ltd. and John Scott
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
 * \file basic-tests.cpp
 * \brief Contains catch testing
 */

// This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <qsl/utils/random.hpp>
#include <qsl/utils/quantum.hpp>
#include <qsl/utils/misc.hpp>
#include <qsl/qubits.hpp>
#include <qsl/verify.hpp>

#include <sstream>

TEST_CASE("Check the measureOut function", "[resize]")
{
    using Sim = qsl::Qubits<qsl::Type::Resize, double>;    

    // Make a simulator in computational basis state
    Sim q{ 4 };
    
    unsigned outcome = q.measureOut(0); 
    REQUIRE(outcome == 0);
    REQUIRE(q.getNumQubits() == 3);

    std::vector<qsl::complex<double>> state = q.getState();
    REQUIRE(state.size() == (1 << 3));
    REQUIRE(qsl::norm(state) - 1 < 1e-10); // Check correct norm
    
    // Check that measuring out a Bell pair works
    Sim q1{ 2 };
    const double sqrt2 = std::sqrt(2);
    q1.setState({{1/sqrt2,0}, {0,0}, {0,0}, {1/sqrt2,0}}); // Set Bell pair
    outcome = q1.measureOut(1); // Measure out qubit 1
    state = q1.getState();
    REQUIRE(state.size() == 2);

    std::vector<qsl::complex<double>> state_correct;
    if (outcome == 0) {
	state_correct.push_back({1,0});
	state_correct.push_back({0,0});
    } else {
	state_correct.push_back({0,0});
	state_correct.push_back({1,0});
    } 
    REQUIRE(distance(state_correct, state) < 1e-10);

    // Check that you can add a qubit
    q1.addQubit();
    REQUIRE(q1.getNumQubits() == 2);
    state = q1.getState();
    // The correct state has new zeroes at the end
    state_correct.push_back({0,0});
    state_correct.push_back({0,0});

    REQUIRE(distance(state_correct, state) < 1e-10);
    REQUIRE(qsl::norm(state) - 1 < 1e-10); // Check correct norm

    
}

TEST_CASE("Test choose function", "[misc-utils]")
{
    REQUIRE(qsl::choose(1,1) == 1);
    REQUIRE(qsl::choose(4,3) == 4);
    REQUIRE(qsl::choose(10,3) == 120);
    REQUIRE(qsl::choose(10,10) == 1);
}

TEST_CASE ("Test overloaded vector substraction function", "[misc-utils]")
{
    const unsigned num_qubits{ 5 };

    const std::vector<qsl::complex<double>> state1
	= qsl::makeRandomState<double>(num_qubits);

    const std::vector<qsl::complex<double>> state2
	= qsl::makeRandomState<double>(num_qubits + 1);

    // Cannot subtract vectors of different sizes
    CHECK_THROWS(state1 - state2);
    
}

TEST_CASE ("Test the convertState function", "[misc-utils]")
{
    const unsigned num_qubits{ 5 };
    
    // Make a random state of doubles
    const std::vector<qsl::complex<double>> state_double
	= qsl::makeRandomState<double>(num_qubits);

    // Convert it to floats
    const std::vector<qsl::complex<float>> state_float
	= qsl::convertState<float>(state_double);

    // Check the two states are equal by verifying all 
    double val = 0;
    for (std::size_t n = 0; n < state_float.size(); n++) {
	val += std::abs(state_double[n].real - state_float[n].real);
	val += std::abs(state_double[n].imag - state_float[n].imag);
    }

    // What is a reasonable number to put here?
    REQUIRE(val < 1e-5);
}

TEST_CASE ("Test the innerProduct function", "[quantum-utils]")
{
    const std::vector<qsl::complex<double>> a{ {1,1.2}, {1.2,1}, {1,1.5} };
    const std::vector<qsl::complex<double>> b{ {1,0}, {0,-1}, {-1,0} };
    const std::vector<qsl::complex<double>> c{
	{1.2,0}, {1.2,-2.1}, {1.3,2}, {0,0} };

    // Check that inner product throws exception for different sized vectors
    CHECK_THROWS(innerProduct(a,c));	 
		 
}


TEST_CASE ("Test the checkStateSize function", "[quantum-utils]")
{
    // Not a valid state size, should throw exception
    const std::vector<qsl::complex<double>> pretend_state{ {1,0}, {0,1}, {1,1} };
    CHECK_THROWS(checkStateSize(pretend_state));

    // Check that the function returns the correct state size
    const std::vector<qsl::complex<double>> state{
	{1,0}, {0,1}, {1,0}, {0,1} };

    REQUIRE( checkStateSize(state) == 2 );
}

TEST_CASE( "Test Fubini Study distance", "[quantum-utils]" )
{
    std::vector<qsl::complex<double>> a{ {1,0}, {0,1}, {1,1} }; 
    std::vector<qsl::complex<double>> double_a{ {2,0}, {0,2}, {2,2} }; 

    std::vector<qsl::complex<double>> b{ {1.3,2}, {1,1.3}, {1.5,0.2} }; 
    std::vector<qsl::complex<double>> double_b{ {2.6,4}, {2,2.6}, {3,0.4} }; 

    // Test that distance between equal and scaled vectors is zero
    REQUIRE( fubiniStudy(a,a) < 1e-13);
    REQUIRE( fubiniStudy(double_a,a) < 1e-13);

    double distance_a_b = fubiniStudy(a,b);
    double scaled_distance_a_b = fubiniStudy(a,double_b);

    REQUIRE( std::abs(distance_a_b - scaled_distance_a_b) < 1e-13);
    
}

TEST_CASE( "Qubits<Default> basic functions", "[state-functions]" )
{
    using Sim = qsl::Qubits<qsl::Type::Default, double>;
    
    const unsigned num_qubits{ 3 };
    Sim q{ num_qubits };
    
    // Set q to a random state
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomState<double>(num_qubits + 1);

    // Can't assign state vector with a different number of qubits
    CHECK_THROWS( q.setState(state) );

    // Check that the basis states work
    std::size_t index = 4; // Should choose a random value
    q.setBasisState(index);
    std::vector<qsl::complex<double>> basis_state = q.getState();

    // Check that the right element is 1
    REQUIRE(basis_state[index].real == 1);
    REQUIRE(std::abs(norm(basis_state) - 1) < 1e-15);
}

TEST_CASE( "Qubits<Omp> basic functions", "[state-functions]" )
{
    using Sim = qsl::Qubits<qsl::Type::Omp, double>;
    
    const unsigned num_qubits{ 3 };
    Sim q{ num_qubits };
    
    // Set q to a random state
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomState<double>(num_qubits + 1);

    // Can't assign state vector with a different number of qubits
    CHECK_THROWS( q.setState(state) );

    // Check that the basis states work
    std::size_t index = 4; // Should choose a random value
    q.setBasisState(index);
    std::vector<qsl::complex<double>> basis_state = q.getState();

    // Check that the right element is 1
    REQUIRE(basis_state[index].real == 1);
    REQUIRE(std::abs(norm(basis_state) - 1) < 1e-15);
}

TEST_CASE( "Qubits<Resize> basic functions", "[state-functions]" )
{
    using Sim = qsl::Qubits<qsl::Type::Resize, double>;
    
    const unsigned num_qubits{ 3 };
    Sim q{ num_qubits };
    
    // Set q to a random state
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomState<double>(num_qubits + 1);

    // Can't assign state vector with a different number of qubits
    CHECK_THROWS( q.setState(state) );

    // Check that the basis states work
    std::size_t index = 4; // Should choose a random value
    q.setBasisState(index);
    std::vector<qsl::complex<double>> basis_state = q.getState();

    // Check that the right element is 1
    REQUIRE(basis_state[index].real == 1);
    REQUIRE(std::abs(norm(basis_state) - 1) < 1e-15);
}


TEST_CASE( "Qubits<NP> basic functions", "[state-functions]" )
{
    using Sim = qsl::Qubits<qsl::Type::NP, double>;

    // Try to make an object with more ones than qubits
    CHECK_THROWS(Sim(4,5));
    
    const unsigned num_qubits{ 4 };
    const unsigned num_ones{ 2 };
    Sim q{ num_qubits };

    // Can't set number of ones more than number of qubits
    CHECK_THROWS(q.setNumOnes(num_qubits + 1));

    // Check that the returned number of ones is the same as that set
    q.setNumOnes(num_ones);
    REQUIRE(q.getNumOnes() == num_ones);
    
    // Set q to a random state
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomNPState<double>(num_qubits + 1);

    // Can't assign state vector with a different number of qubits
    CHECK_THROWS( q.setState(state) );

    // Check that the basis states work
    std::size_t index = 4; // Should choose a random value
    q.setBasisState(index);
    std::vector<qsl::complex<double>> basis_state = q.getState();

    // Check that the right element is 1
    REQUIRE(basis_state[index].real == 1);
    REQUIRE(std::abs(norm(basis_state) - 1) < 1e-15);
}



TEST_CASE( "Qubits<Default> object constructors", "[constructors]" )
{
    using Sim = qsl::Qubits<qsl::Type::Default, double>;

    const unsigned num_qubits{ 5 };
    Sim q{ num_qubits };

    // Set q to a random state
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomState<double>(num_qubits);
    q.setState(state);

    // Construct another object using the copy constructor
    Sim q_copy{ q };    

    // Check that the other properties are also the same
    REQUIRE(q.getNumQubits() == q_copy.getNumQubits());
    
    // Check that both objects have the same internal state
    double distance = qsl::distance(q.getState(), q_copy.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Construct another object using initialisation from state vector
    Sim q_from_state{ state };    

    // Check that both objects have the same internal state
    distance = qsl::distance(q.getState(), q_from_state.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Check that assignment works
    Sim q_new{ num_qubits };
    q_new = q; // Perform the test
    distance = qsl::distance(q.getState(), q_new.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Check that assignment fails for wrong number of qubits
    Sim q_wrong{ num_qubits + 1};
    CHECK_THROWS(q_wrong = q); // Perform the test

}

TEST_CASE( "Qubits<OMP> object constructors", "[constructors]" )
{
    using Sim = qsl::Qubits<qsl::Type::Omp, double>;

    const unsigned num_qubits{ 5 };
    Sim q{ num_qubits };

    // Set q to a random state
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomState<double>(num_qubits);
    q.setState(state);

    // Construct another object using the copy constructor
    Sim q_copy{ q };    

    // Check that the other properties are also the same
    REQUIRE(q.getNumQubits() == q_copy.getNumQubits());
    
    // Check that both objects have the same internal state
    double distance = qsl::distance(q.getState(), q_copy.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Construct another object using initialisation from state vector
    Sim q_from_state{ state };    

    // Check that both objects have the same internal state
    distance = qsl::distance(q.getState(), q_from_state.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Check that assignment works
    Sim q_new{ num_qubits };
    q_new = q; // Perform the test
    distance = qsl::distance(q.getState(), q_new.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Check that assignment fails for wrong number of qubits
    Sim q_wrong{ num_qubits + 1};
    CHECK_THROWS(q_wrong = q); // Perform the test
}

TEST_CASE( "Qubits<NP> object constructors", "[constructors]" )
{
    using Sim = qsl::Qubits<qsl::Type::NP, double>;

    const unsigned num_qubits{ 5 };
    Sim q{ num_qubits };

    // Set q to a random state (number preserving
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomNPState<double>(num_qubits);
    q.setState(state);

    // Construct another object using the copy constructor
    Sim q_copy{ q };    

    // Check that the other properties are also the same
    REQUIRE(q.getNumQubits() == q_copy.getNumQubits());
    
    // Check that both objects have the same internal state
    double distance = qsl::distance(q.getState(), q_copy.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Construct another object using initialisation from state vector
    Sim q_from_state{ state };    

    // Check that both objects have the same internal state
    distance = qsl::distance(q.getState(), q_from_state.getState());
    REQUIRE(std::abs(distance) < 1e-10);
    
    // Check that assignment works
    Sim q_new{ num_qubits };
    q_new = q; // Perform the test
    distance = qsl::distance(q.getState(), q_new.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Check that assignment fails for wrong number of qubits
    Sim q_wrong{ num_qubits + 1};
    CHECK_THROWS(q_wrong = q); // Perform the test

}

/// NOT TESTING ANYTHING YET
TEST_CASE( "Qubits<Default> against Qubits<OMP>", "[compare]" )
{
    using Sim1 = qsl::Qubits<qsl::Type::Default, double>;
    using Sim2 = qsl::Qubits<qsl::Type::Omp, double>;

    qsl::Verify<Sim1, Sim2, qsl::DefaultStateGen<double>,
		qsl::MeasureChecker,
		qsl::PostselectChecker,
		qsl::SampleChecker,
		qsl::SampleAllChecker,
		qsl::ProbChecker,
		qsl::DefaultGateChecker> verify;
    verify.configureState(8);
    verify.configureChecker<qsl::SampleAllChecker>(1000000, 0.99);
    //verify.checkAll();

    auto prob_result = verify.check<qsl::ProbChecker>();
    // Check that the probabilities are equal
    REQUIRE(std::abs(prob_result.mean<double>(1)) < 1e-10);
    REQUIRE(prob_result.variance<double>(1) < 1e-10);
    REQUIRE(std::abs(prob_result.mean<double>(2)) < 1e-10);
    REQUIRE(prob_result.variance<double>(2) < 1e-10);
    
    // Check the measurement
    auto measure_result = verify.check<qsl::MeasureChecker>();
    // Check that all the intervals overlap
    REQUIRE(std::abs(measure_result.mean<double>(1) - 1) < 1e-10);
    REQUIRE(measure_result.variance<double>(1) < 1e-10);

    auto postselect_result = verify.check<qsl::PostselectChecker>();
    // Check that all the probability differences are zero for both outcomes
    REQUIRE(std::abs(postselect_result.mean<double>(1)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(1) < 1e-10);
    REQUIRE(std::abs(postselect_result.mean<double>(2)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(2) < 1e-10);
    // Check that all the state distances are zero for both outcomes
    REQUIRE(std::abs(postselect_result.mean<double>(3)) < 1e-7);
    REQUIRE(postselect_result.variance<double>(3) < 1e-7);
    REQUIRE(std::abs(postselect_result.mean<double>(4)) < 1e-7);
    REQUIRE(postselect_result.variance<double>(4) < 1e-7);
    
    auto sample_result = verify.check<qsl::SampleChecker>();
    // Check that all the confidence intervals overlap
    REQUIRE(std::abs(sample_result.mean<double>(1) - 1) < 1e-10);
    
    auto sample_all_result = verify.check<qsl::SampleAllChecker>();
    // Check that all the confidence intervals overlap
    REQUIRE(std::abs(sample_all_result.overlap_rate - 1) < 1e-10);
    // Check the mean relative error is less than 15%
    REQUIRE(std::abs(sample_all_result.mean_relative_error) < 0.15);
        
    auto gate_result = verify.check<qsl::DefaultGateChecker>();
    for (const auto & [gate,table] : gate_result) {
	// Check that all the distances are zero
	REQUIRE(std::abs(table.mean<double>("distance")) < 1e-9);
	REQUIRE(std::abs(table.variance<double>("distance")) < 1e-9);
    }
}

TEST_CASE( "Qubits<Default> against Qubits<NP>", "[compare]" )
{
    using Sim1 = qsl::Qubits<qsl::Type::Default, double>;
    using Sim2 = qsl::Qubits<qsl::Type::NP, double>;

    qsl::Verify<Sim1, Sim2, qsl::NPStateGen<double>,
		qsl::MeasureChecker,
		qsl::PostselectChecker,
		qsl::SampleChecker,
		qsl::SampleAllChecker,
		qsl::ProbChecker,
		qsl::NPGateChecker> verify;
    verify.configureState(8, 5); // Pass number of ones as second argument
    verify.configureChecker<qsl::SampleAllChecker>(1000000, 0.99);
    //verify.checkAll();

    // Check the measurement
    auto measure_result = verify.check<qsl::MeasureChecker>();
    // Check that all the intervals overlap
    REQUIRE(std::abs(measure_result.mean<double>(1) - 1) < 1e-10);
    REQUIRE(measure_result.variance<double>(1) < 1e-10);

    auto postselect_result = verify.check<qsl::PostselectChecker>();
    // Check that all the probability differences are zero for both outcomes
    REQUIRE(std::abs(postselect_result.mean<double>(1)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(1) < 1e-10);
    REQUIRE(std::abs(postselect_result.mean<double>(2)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(2) < 1e-10);
    // Check that all the state distances are zero for both outcomes
    REQUIRE(std::abs(postselect_result.mean<double>(3)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(3) < 1e-10);
    REQUIRE(std::abs(postselect_result.mean<double>(4)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(4) < 1e-10);
    
    auto sample_result = verify.check<qsl::SampleChecker>();
    // Check that all the confidence intervals overlap
    REQUIRE(std::abs(sample_result.mean<double>(1) - 1) < 1e-10);
    
    auto sample_all_result = verify.check<qsl::SampleAllChecker>();
    // Check that all the confidence intervals overlap
    REQUIRE(std::abs(sample_all_result.overlap_rate - 1) < 1e-10);
    // Check the mean relative error is less than 15%
    REQUIRE(std::abs(sample_all_result.mean_relative_error) < 0.15);
    
    auto prob_result = verify.check<qsl::ProbChecker>();
    // Check that the probabilities are equal
    REQUIRE(std::abs(prob_result.mean<double>(1)) < 1e-10);
    REQUIRE(prob_result.variance<double>(1) < 1e-10);
    REQUIRE(std::abs(prob_result.mean<double>(2)) < 1e-10);
    REQUIRE(prob_result.variance<double>(2) < 1e-10);
    
    auto gate_result = verify.check<qsl::NPGateChecker>();
    for (const auto & [gate,table] : gate_result) {
	// Check that all the distances are zero
	REQUIRE(std::abs(table.mean<double>("distance")) < 1e-9);
	REQUIRE(std::abs(table.variance<double>("distance")) < 1e-9);
    }
}

TEST_CASE( "Qubits<Default> against Qubits<Resize>", "[compare]" )
{
    using Sim1 = qsl::Qubits<qsl::Type::Default, double>;
    using Sim2 = qsl::Qubits<qsl::Type::Resize, double>;

    qsl::Verify<Sim1, Sim2, qsl::DefaultStateGen<double>,
		qsl::MeasureChecker,
		qsl::PostselectChecker,
		qsl::SampleChecker,
		qsl::SampleAllChecker,
		qsl::ProbChecker,
		qsl::DefaultGateChecker> verify;
    verify.configureState(8);
    verify.configureChecker<qsl::SampleAllChecker>(1000000, 0.99);
    //verify.checkAll();

    // Check the measurement
    auto measure_result = verify.check<qsl::MeasureChecker>();
    // Check that all the intervals overlap
    REQUIRE(std::abs(measure_result.mean<double>(1) - 1) < 1e-10);
    REQUIRE(measure_result.variance<double>(1) < 1e-10);

    auto postselect_result = verify.check<qsl::PostselectChecker>();
    // Check that all the probability differences are zero for both outcomes
    REQUIRE(std::abs(postselect_result.mean<double>(1)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(1) < 1e-10);
    REQUIRE(std::abs(postselect_result.mean<double>(2)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(2) < 1e-10);
    // Check that all the state distances are zero for both outcomes
    REQUIRE(std::abs(postselect_result.mean<double>(3)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(3) < 1e-10);
    REQUIRE(std::abs(postselect_result.mean<double>(4)) < 1e-10);
    REQUIRE(postselect_result.variance<double>(4) < 1e-10);
    
    auto sample_result = verify.check<qsl::SampleChecker>();
    // Check that all the confidence intervals overlap
    REQUIRE(std::abs(sample_result.mean<double>(1) - 1) < 1e-10);
    
    auto sample_all_result = verify.check<qsl::SampleAllChecker>();
    // Check that all the confidence intervals overlap
    REQUIRE(std::abs(sample_all_result.overlap_rate - 1) < 1e-10);
    // Check the mean relative error is less than 15%
    REQUIRE(std::abs(sample_all_result.mean_relative_error) < 0.15);
    
    auto prob_result = verify.check<qsl::ProbChecker>();
    // Check that the probabilities are equal
    REQUIRE(std::abs(prob_result.mean<double>(1)) < 1e-10);
    REQUIRE(prob_result.variance<double>(1) < 1e-10);
    REQUIRE(std::abs(prob_result.mean<double>(2)) < 1e-10);
    REQUIRE(prob_result.variance<double>(2) < 1e-10);
    
    auto gate_result = verify.check<qsl::DefaultGateChecker>();
    for (const auto & [gate,table] : gate_result) {
	// Check that all the distances are zero
	REQUIRE(std::abs(table.mean<double>("distance")) < 1e-9);
	REQUIRE(std::abs(table.variance<double>("distance")) < 1e-9);
    }
}
