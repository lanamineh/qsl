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
