/**
 * \file qsl-testing.cpp
 *
 * Contains gtest code.
 */

#include <gtest/gtest.h>

#include "meta.hpp"
#include "qsl.hpp"
#include <armadillo>

// Make combination of all simulator arguments
using sim_args = qsl::comb_t<std::tuple<float, double, long double>,
			     std::tuple<std::true_type, std::false_type>,
			     std::tuple<qsl::seq, qsl::omp, qsl::opt>>;
using test_sims = qsl::gtest_typelist_t<qsl::cat_t<qsl::make_sims_t<qsl::basic, sim_args>,
						   qsl::make_sims_t<qsl::resize, sim_args>,
						   qsl::make_sims_t<qsl::number, sim_args>>>;


template <typename T>
class ConstructorTests : public testing::Test {};

TYPED_TEST_SUITE(ConstructorTests, test_sims);

/**
 * \brief Testing the default constuctor.
 *
 * Creates an empty simulator, checks the number of qubits is zero, that the
 * size of the state vector is 1 and that it is normalised. It tests both operator[]
 * and the .state() method.
 */
TYPED_TEST(ConstructorTests, DefaultConstructor)
{
    using sim_t = TypeParam;
    
    sim_t q;
    EXPECT_EQ(q.qubits(), 0) << "There must be 0 qubits.";
    EXPECT_EQ(q.size(), 1) << "The state vector size must be 1.";
    EXPECT_NEAR(std::abs(q[0]), 1, 1e-7) << "The state vector must be normalised.";
}

/**
 * \brief Test the qubit number constructor.
 * 
 * Try to initialise a state vector of 0 to 20 qubits, each time checking
 * qubits number, state vector size and initialised to all zeros. For debug mode,
 * try to create a very large state vector.
 */
TYPED_TEST(ConstructorTests, QubitConstructor)
{
    using sim_t = TypeParam;
    
    for (unsigned i = 0; i < 21; i++) {
	sim_t q{i};
	EXPECT_EQ(q.qubits(), i) << "There must be " << i << " qubits.";
	EXPECT_EQ(q.size(), 1 << i) << "The state vector size must be " << (1 << i);
	EXPECT_NEAR(std::abs(q[0]), 1, 1e-7) << "The state vector must be normalised.";
	for (std::size_t j = 1; j < q.size(); j++) {
	    EXPECT_NEAR(std::abs(q[j]), 0, 1e-7) << "The state vector must be normalised.";
	}
    }

    if constexpr (qsl::is_debug<sim_t>)
    {
	EXPECT_THROW(sim_t q{100}, std::runtime_error) << "Should not be able to allocate 100 qubits.";
    }
}
