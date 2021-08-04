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
 * \file random.hpp
 * \brief A class for generating random numbers, states, etc.
 */

#ifndef QSL_RANDOM_HPP
#define QSL_RANDOM_HPP


#include <random>
#include "complex.hpp"

namespace qsl {

    /**
     * \brief A small class for generating random numbers
     *
     * This class is a wrapper around a pseudo-random 
     * number generator. Once the class has been constructed,
     * use the getNum() member function to get a uniform
     * random number in a range specified in the constructor.
     * 
     * The implementation uses the Mersenne-Twister algorithm
     * with n=19937 (called mt19937), which passes enough 
     * statistical tests for randomness that it is legitimate
     * to use it in scientific software.
     *
     * std::default_random_engine is not a good choice because
     * it is designed for "relatively casual, inexpert and/or
     * lightweight use" (C++ std library documentation).
     *
     */
    template<std::floating_point Fp = double>
    class Random
    {
	std::mt19937 generator;
	std::uniform_real_distribution<Fp> distribution;

    public:

	/**
	 * \brief Construct the class and specify the range
	 *
	 * Constructs a Random object which can be used
	 * to get uniform random numbers in the range
	 * [a,b).
	 *
	 * \param a The bottom end of the range
	 * \param b The top end of the range
	 */
	Random(Fp a, Fp b);
    
	/**
	 * \brief Get a uniform random number in the range [a,b)
	 *
	 * The range [a,b) is specified in the constructor.
	 *
	 */
	Fp getNum();
    };

    /// Deduction guide for Random class
    //template<std::floating_point Fp = double>
    //Random(Fp, Fp) -> Random<Fp>;



    /**
     * \brief Restrictions that apply to some simulators
     *
     * This class acts as a tag that distinguishes between
     * different simulators which have different capabilities.
     * None means that there are no restrictions. NumberPreserving
     * simulators are only allowed a fixed number of ones.
     */
    enum class Restrictions
    {
	None,
	NumberPreserved,
    };

    /**
     * \brief Make a random state vector
     */
    template<std::floating_point Fp = double>
    std::vector<complex<Fp>> makeRandomState(unsigned nqubits);

    /**
     * \brief Make a random number preserved state vector
     *
     * Makes a random number preserved state, with a random 
     * number of ones between 1 and (nqubits-1) inclusive.
     * 
     */
    template<std::floating_point Fp = double>
    std::vector<complex<Fp>> makeRandomNPState(unsigned nqubits);

    template<std::floating_point Fp = double>
    std::vector<complex<Fp>> makeRandomNPState(unsigned nqubits, unsigned nones);

    /**
     * \brief Generate a random state vector
     *
     *
     */
    template<std::floating_point Fp, Restrictions R> class RandomStateGen;

    /// No restrictions
    template<std::floating_point Fp>
    class RandomStateGen<Fp, Restrictions::None>
    {
	const std::vector<complex<Fp>> state;
    public:    
	std::vector<complex<Fp>> get() const { return state; };
	RandomStateGen(unsigned nqubits)
	    : state{ makeRandomState<Fp>(nqubits) }
	    { }
    };

    /// Number preserved
    template<std::floating_point Fp>
    class RandomStateGen<Fp, Restrictions::NumberPreserved>
    {
	const std::vector<complex<Fp>> state;
    public:
	std::vector<complex<Fp>> get() const { return state; };
	RandomStateGen(unsigned nqubits)
	    : state{ makeRandomNPState<Fp>(nqubits) }
	    { }
    };

    /**
     * \brief Make a vector of random phases (doubles)
     *
     * This is used in the benchmarking and verification
     * tests
     *
     */
    template<std::floating_point Fp = double>
    std::vector<Fp> makeRandomPhases(int test_len);

}
#endif
