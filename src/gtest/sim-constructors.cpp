#include <gtest/gtest.h>

#include <qsl/qubits.hpp>
#include <qsl/utils.hpp>
#include <list>

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
class SimBasics : public testing::Test {
public:
    using List = std::list<T>;
    static T shared_;
    T value_;
};

template <typename T>
class NPSimBasics : public testing::Test {
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



/**
 * \brief Declare the test suite that depends on these types
 *
 * To use the test suite, write TYPED_TEST(FpUtilities...
 * and then use TypeParam wherever you would use a type from the
 * FpTypes list. The test suite will automatically be performed
 * for every type in the list.
 *
 */
TYPED_TEST_SUITE(SimBasics, SimTypes);

TYPED_TEST_SUITE(NPSimBasics, NPSimTypes);


/// Qubits object constructors
TYPED_TEST(SimBasics, ConstructorTest)
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

    // Check that assignment fails for wrong number of qubits
    Sim q_wrong{ num_qubits + 1};
    EXPECT_THROW(q_wrong = q, std::logic_error); // Perform the test
}

// Test for NP constructors
TYPED_TEST(NPSimBasics, ConstructorTest)
{
    using Sim = TypeParam::Sim;
    using Fp_type = TypeParam::Sim::Fp_type;

    const unsigned num_qubits{ 5 };
    Sim q{ num_qubits };

    // Set q to a random state (number preserving
    const std::vector<qsl::complex<Fp_type>> state
	= qsl::makeRandomNPState<Fp_type>(num_qubits);
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

    // Check that assignment fails for wrong number of qubits
    Sim q_wrong{ num_qubits + 1};
    EXPECT_THROW(q_wrong = q, std::logic_error); // Perform the test

}

/// Test the basic utility functions of the simulators
TYPED_TEST(SimBasics, UtilityFunctions)
{
    using Sim = TypeParam::Sim;
    using Fp_type = TypeParam::Sim::Fp_type;
    
    const unsigned num_qubits{ 3 };
    Sim q{ num_qubits };
    
    // Set q to a random state
    const std::vector<qsl::complex<Fp_type>> state
	= qsl::makeRandomState<Fp_type>(num_qubits + 1);

    // Can't assign state vector with a different number of qubits
    EXPECT_THROW( q.setState(state), std::logic_error);

    // Check that the basis states work
    std::size_t index = 4; // Should choose a random value
    q.setBasisState(index);
    std::vector<qsl::complex<Fp_type>> basis_state = q.getState();

    // Check that the right element is 1
    EXPECT_FLOAT_EQ(basis_state[index].real, 1);
    EXPECT_FLOAT_EQ(std::abs(norm(basis_state)), 1);
}

/// Test number preserving simulator utility functions
TYPED_TEST(NPSimBasics, UtilityFunctions)
{
    using Sim = TypeParam::Sim;
    using Fp_type = TypeParam::Sim::Fp_type;

    // Try to make an object with more ones than qubits
    EXPECT_THROW(Sim(4,5), std::logic_error);
    
    const unsigned num_qubits{ 4 };
    const unsigned num_ones{ 2 };
    Sim q{ num_qubits };

    // Can't set number of ones more than number of qubits
    EXPECT_THROW(q.setNumOnes(num_qubits + 1), std::logic_error);

    // Check that the returned number of ones is the same as that set
    q.setNumOnes(num_ones);
    EXPECT_EQ(q.getNumOnes(), num_ones);
    
    // Set q to a random state
    const std::vector<qsl::complex<Fp_type>> state
	= qsl::makeRandomNPState<Fp_type>(num_qubits + 1);

    // Can't assign state vector with a different number of qubits
    EXPECT_THROW(q.setState(state), std::logic_error);

    // Check that the basis states work
    std::size_t index = 4; // Should choose a random value
    q.setBasisState(index);
    std::vector<qsl::complex<Fp_type>> basis_state = q.getState();

    // Check that the right element is 1
    EXPECT_FLOAT_EQ(basis_state[index].real, 1);
    EXPECT_FLOAT_EQ(std::abs(norm(basis_state)), 1);
}

/// Check the measure out function of the resize qubit
TEST(ResizeSim, MeasureOutTest)
{
    using Sim = qsl::Qubits<qsl::Type::Resize, double>;    

    // Make a simulator in computational basis state
    Sim q{ 4 };
    
    unsigned outcome = q.measureOut(0); 
    EXPECT_EQ(outcome, 0);
    EXPECT_EQ(q.getNumQubits(), 3);

    std::vector<qsl::complex<double>> state = q.getState();
    EXPECT_EQ(state.size(), (1 << 3));
    EXPECT_FLOAT_EQ(qsl::norm(state), 1); // Check correct norm
    
    // Check that measuring out a Bell pair works
    Sim q1{ 2 };
    const double sqrt2 = std::sqrt(2);
    q1.setState({{1/sqrt2,0}, {0,0}, {0,0}, {1/sqrt2,0}}); // Set Bell pair
    outcome = q1.measureOut(1); // Measure out qubit 1
    state = q1.getState();
    EXPECT_EQ(state.size(), 2);

    std::vector<qsl::complex<double>> state_correct;
    if (outcome == 0) {
	state_correct.push_back({1,0});
	state_correct.push_back({0,0});
    } else {
	state_correct.push_back({0,0});
	state_correct.push_back({1,0});
    } 
    EXPECT_FLOAT_EQ(distance(state_correct, state), 0);

    // Check that you can add a qubit
    q1.addQubit();
    EXPECT_EQ(q1.getNumQubits(), 2);
    state = q1.getState();
    // The correct state has new zeroes at the end
    state_correct.push_back({0,0});
    state_correct.push_back({0,0});

    EXPECT_FLOAT_EQ(distance(state_correct, state), 0);
    EXPECT_FLOAT_EQ(qsl::norm(state), 1); // Check correct norm

    
}
