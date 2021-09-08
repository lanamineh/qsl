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
 * \file utils.cpp
 * \brief Contains tests for utility functions
 */

#include <gtest/gtest.h>

#include <qsl/utils.hpp>
#include <list>
#include <sstream>

/**
 * \brief Typed test suite for the floating point utilities
 *
 * This test suite checks the utilities that depend on a floating
 * point parameter for both float and double.
 */
template <typename T>
class FpUtilities : public testing::Test {
public:
    using List = std::list<T>;
    static T shared_;
    T value_;
};

/// List of floating point types to check
using FpTypes = ::testing::Types<double, float>;

/**
 * \brief Declare the test suite that depends on these types
 *
 * To use the test suite, write TYPED_TEST(FpUtilities...
 * and then use TypeParam wherever you would use a type from the
 * FpTypes list. The test suite will automatically be performed
 * for every type in the list.
 *
 */
TYPED_TEST_SUITE(FpUtilities, FpTypes);

/// Test the combinations function 
TEST(Utilities, ChooseFunctionTest)
{
    EXPECT_EQ(qsl::choose(1,1), 1);
    EXPECT_EQ(qsl::choose(4,3), 4);
    EXPECT_EQ(qsl::choose(10,3), 120);
    EXPECT_EQ(qsl::choose(10,10), 1);
}

/// Test overloaded vector substraction function
TYPED_TEST(FpUtilities, VectorSubtraction)
{
    const unsigned num_qubits{ 5 };

    const std::vector<qsl::complex<TypeParam>> state1
	= qsl::makeRandomState<TypeParam>(num_qubits);

    const std::vector<qsl::complex<TypeParam>> state2
	= qsl::makeRandomState<TypeParam>(num_qubits + 1);

    // Cannot subtract vectors of different sizes
    EXPECT_THROW(state1 - state2, std::logic_error);  
}

/// Test the convertState function
TEST(Utilities, ConvertStateTest)
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
    EXPECT_NEAR(val,0,1e-5);
}

/// Test the convertVector function
TEST(Utilities, ConvertVectorTest)
{
    // Make a test vector
    const std::vector<double> vec{ 1,3.2,3,4.9,3.2 }; 
    const std::vector<float> copy{ qsl::convertVector<float>(vec) };

    // Check that the two are equal
    EXPECT_EQ(copy.size(), vec.size());
    for (std::size_t n = 0; n < vec.size(); n++) {
	EXPECT_FLOAT_EQ(vec[n], copy[n]);
    }
    
}


/// Test the innerProduct function
TYPED_TEST(FpUtilities, InnerProductTest)
{
    const std::vector<qsl::complex<TypeParam>> a{ {1,1.2}, {1.2,1}, {1,1.5} };
    const std::vector<qsl::complex<TypeParam>> b{ {1,0}, {0,-1}, {-1,0} };
    const std::vector<qsl::complex<TypeParam>> c{
	{1.2,0}, {1.2,-2.1}, {1.3,2}, {0,0} };

    // Check that inner product throws exception for different sized vectors
    EXPECT_THROW(innerProduct(a,c), std::logic_error);	 
		 
}

/// Test the checkStateSize function
TYPED_TEST(FpUtilities, CheckStateSizeTest)
{
    // Not a valid state size, should throw exception
    const std::vector<qsl::complex<TypeParam>> pretend_state{ {1,0}, {0,1}, {1,1} };
    EXPECT_THROW(checkStateSize(pretend_state), std::logic_error);

    // Check that the function returns the correct state size
    const std::vector<qsl::complex<TypeParam>> state{
	{1,0}, {0,1}, {1,0}, {0,1} };

    EXPECT_EQ(checkStateSize(state), 2);
}

/// Test Fubini Study distance
TYPED_TEST(FpUtilities, FubiniStudyTest)
{
    std::vector<qsl::complex<TypeParam>> a{ {1,0}, {0,1}, {1,1} }; 
    std::vector<qsl::complex<TypeParam>> TypeParam_a{ {2,0}, {0,2}, {2,2} }; 

    std::vector<qsl::complex<TypeParam>> b{ {1.3,2}, {1,1.3}, {1.5,0.2} }; 
    std::vector<qsl::complex<TypeParam>> TypeParam_b{ {2.6,4}, {2,2.6}, {3,0.4} }; 

    // Test that distance between equal and scaled vectors is zero
    EXPECT_NEAR(fubiniStudy(a,a), 0, 1e-13);
    EXPECT_NEAR(fubiniStudy(TypeParam_a,a), 0, 1e-13);

    TypeParam distance_a_b = fubiniStudy(a,b);
    TypeParam scaled_distance_a_b = fubiniStudy(a,TypeParam_b);

    EXPECT_NEAR(std::abs(distance_a_b - scaled_distance_a_b), 0, 1e-13);
    
}

/// Check the next (number with fixed number of ones) function
TEST(Utilities, NextFunctionTest)
{
    std::size_t x{ 0b111 };
    qsl::next(x);
    EXPECT_EQ(x,0b1011);
    qsl::next(x);
    EXPECT_EQ(x,0b1101);
    qsl::next(x);
    EXPECT_EQ(x,0b1110);
    qsl::next(x);
    EXPECT_EQ(x,0b10011);
    qsl::next(x);
    EXPECT_EQ(x,0b10101);
    qsl::next(x);
    EXPECT_EQ(x,0b10110);
    qsl::next(x);
    EXPECT_EQ(x,0b11001);
    qsl::next(x);
    EXPECT_EQ(x,0b11010);
    qsl::next(x);
    EXPECT_EQ(x,0b11100);
}

/// Test the hamming weight function
TEST(Utilities, HammingWeightTest)
{
    EXPECT_EQ(qsl::hammingWeight(0b11011010), 5);
    EXPECT_EQ(qsl::hammingWeight(0b1010), 2);
    EXPECT_EQ(qsl::hammingWeight(0b1), 1);
    EXPECT_EQ(qsl::hammingWeight(0b0), 0);
}

/// Test the complex struct
TYPED_TEST(FpUtilities, ComplexStructTest)
{
    // Check explicit assignment
    qsl::complex<TypeParam> x0{ .real = 1, .imag = 3};
    EXPECT_EQ(x0.real, 1);
    EXPECT_EQ(x0.imag, 3);

    // Check real,imag constructor
    qsl::complex<TypeParam> x1{-1,4.5};
    EXPECT_EQ(x1.real, -1);
    EXPECT_EQ(x1.imag, 4.5);

    // Check default constructor
    qsl::complex<TypeParam> x2;
    EXPECT_EQ(x2.real, 0);
    EXPECT_EQ(x2.imag, 0);

    // Check substraction
    EXPECT_EQ((x0 - x1).real, 2);
    EXPECT_EQ((x0 - x1).imag, -1.5);
    
}

/// Test absolute value function
TYPED_TEST(FpUtilities, AbsTest)
{
    qsl::complex<TypeParam> x0;
    EXPECT_FLOAT_EQ(qsl::abs(x0), 0);

    qsl::complex<TypeParam> x1{1,0};
    EXPECT_FLOAT_EQ(qsl::abs(x1), 1);

    qsl::complex<TypeParam> x2{0,1};
    EXPECT_FLOAT_EQ(qsl::abs(x2), 1);

    qsl::complex<TypeParam> x3{0,-1};
    EXPECT_FLOAT_EQ(qsl::abs(x3), 1);

    qsl::complex<TypeParam> x4{1,-1};
    EXPECT_FLOAT_EQ(qsl::abs(x4), std::sqrt(2));
  
}

/// Test the random number generator class
TYPED_TEST(FpUtilities, RandomClassTest)
{
    // Check that a zero length range works
    qsl::Random<TypeParam> rand{1.2,1.2};

    EXPECT_FLOAT_EQ(rand.getNum(), 1.2);

    ///\todo How to check that the distribution is uniform and random?
    
}

/// Test the random state generator class
TYPED_TEST(FpUtilities, MakeRandomStateTest)
{
    // Number of qubits
    const unsigned nqubits{ 5 };
    
    // Check that a zero length range works
    std::vector<qsl::complex<TypeParam>> state{
	qsl::makeRandomState<TypeParam>(nqubits)
    };

    // Check there are the correct number of elements
    EXPECT_EQ(state.size(), 1 << nqubits);

    // Check that the state is normalised
    EXPECT_FLOAT_EQ(qsl::norm(state), 1);
    
}

/// Test the function which test for number-preserving states
TYPED_TEST(FpUtilities, CheckStateNPTest)
{
    // Number-preserved state
    std::vector<qsl::complex<TypeParam>> s0 {
	{0,0}, {0,1}, {1,0}, {0,0}
    };
    qsl::normalise(s0);
    EXPECT_EQ(checkStateNP(s0), 1);

    // Edge-case no ones
    std::vector<qsl::complex<TypeParam>> s1 {
	{1,0}, {0,0}, {0,0}, {0,0}
    };
    qsl::normalise(s1);
    EXPECT_EQ(checkStateNP(s1), 0);

    // Edge-case all ones
    std::vector<qsl::complex<TypeParam>> s2 {
	{0,0}, {0,0}, {0,0}, {1,0}
    };
    qsl::normalise(s2);
    EXPECT_EQ(checkStateNP(s2), 2);
    
    // Non-number-preserved state
    std::vector<qsl::complex<TypeParam>> s3 {
	{1,0}, {0,1}, {1,0}, {0,0}
    };
    qsl::normalise(s3);
    EXPECT_THROW(checkStateNP(s3), std::logic_error);

}

/// Test the random number-preserving state generator function
TYPED_TEST(FpUtilities, MakeRandomNPStateTest)
{
    // Number of qubits
    const unsigned nqubits{ 5 };
    const unsigned nones{ 3 };
    
    // Make a random number preserved state
    std::vector<qsl::complex<TypeParam>> state{
	qsl::makeRandomNPState<TypeParam>(nqubits, nones)
    };
    EXPECT_EQ(state.size(), 1 << nqubits); // Correct dimension
    EXPECT_FLOAT_EQ(qsl::norm(state), 1); // Is normalised
    EXPECT_EQ(checkStateNP(state), nones); // Correct number of ones

    // Check that random number of ones variant
    std::vector<qsl::complex<TypeParam>> state_random{
	qsl::makeRandomNPState<TypeParam>(nqubits)
    };
    EXPECT_EQ(state_random.size(), 1 << nqubits); // Correct dimension
    EXPECT_FLOAT_EQ(qsl::norm(state_random), 1); // Is normalised
    EXPECT_LE(checkStateNP(state_random), nqubits); // Is NP and nones <= nqubits
    
}

/// Check the random phases function
TYPED_TEST(FpUtilities, RandomPhasesTest)
{
    const std::size_t length{ 35 };
    
    // Make a list of random phases
    std::vector<TypeParam> phases{
	qsl::makeRandomPhases<TypeParam>(length)
    };
    EXPECT_EQ(phases.size(), length); // Check correct length
    for (std::size_t n = 0; n < length; n++) {
	// Check elements are within the range -pi to pi
	EXPECT_LE(phases[n], M_PI);
	EXPECT_GE(phases[n], -M_PI);
    }


}

/// Test overloaded vector printing
TEST(Utilities, VectorPrinting)
{
    std::stringstream ss1, ss2;
    std::vector<int> vec { 1,2,5,4,4 };
    ss1 << vec;
    // Check the correct thing is printed
    ss2 << 1 << std::endl
	<< 2 << std::endl
	<< 5 << std::endl
	<< 4 << std::endl
	<< 4 << std::endl;

    EXPECT_EQ(ss1.str(), ss2.str());
}
