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
 * \file misc.hpp
 * \brief Contains miscellaneous useful functions (to be sorted)
 *
 */

#ifndef QSL_MISC_HPP
#define QSL_MISC_HPP

#include <ostream>
#include <vector>
#include <random>
#include "complex.hpp"

namespace qsl {

    /**
     *
     * \brief Get the next number with the same number of bits set.
     *
     * This function will get the next number with the same number
     * of ones as x. The result is stored in x.
     *
     * This is from here:
     * https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
     * 
     */
    void next(std::size_t & x);

    /**
     * \brief Print the timing results for a benchmark
     */
    void printRow(std::string col1, double time1, double time2);

    /**
     * \brief Print the timing results for a benchmark (no column name)
     */
    void printRow(double time1, double time2);

    /**
     * \brief Print a number in binary with a label
     *
     * \todo The length of the bitstring is hardcoded, fix it
     */
    void printBin(const std::string & name, unsigned int x);

    /**
     * \brief Compute the difference between two vector
     *
     * This function is not meant to be fast.
     */
    template<typename T>
    std::vector<T> operator - (const std::vector<T> & v, const std::vector<T> & w)
    {    
	if(v.size() != w.size()) {
	    std::string msg = "Cannot subtract vectors of different sizes: ";
	    msg += std::to_string(v.size()) + " and " + std::to_string(w.size());
	    throw std::logic_error(msg);
	}

	std::vector<T> result;
	for(std::size_t n=0; n<v.size(); n++) {
	    result.push_back(v[n] - w[n]);
	}

	return result;   
    }

    /**
     * \brief Compute the distance between two real or complex vectors
     *
     * This function computes the distance between two vectors,
     * which is the norm of their difference. This function is not
     * meant to be fast
     *
     */
    template<typename T>
    double distance(const std::vector<T> & v1, const std::vector<T> & v2)
    {
	return norm(v1 - v2);
    }


    /**
     * \brief Convert a state vector from one floating point type to another
     *
     * This function converts a state vector whose amplitudes use one floating
     * point precision (such as double) to a different floating point precision
     * (such as float). The function is not useful in the normal use of any of
     * the simulators, but is necessary when benchmarking or verifying simulators
     * with different precisions. For example, it might be necessary to 
     * instantiate two simulators from the same random state vector with
     * double precision. If one of the simulators uses floats, then the random
     * state vector will need to be converted to float.
     *
     * \param state The input state vector to be converted
     * \return The state vector in the new floating point precision
     *
     * \todo At the moment, template type checking is performed using SFINAE
     * based on std::enable_if (checking if From is convertible to To). Find
     * a more modern alternative
     *
     */
    template<typename To, typename From>
    typename std::enable_if_t<std::is_convertible_v<From, To>,
			      std::vector<complex<To>>>
    convertState(const std::vector<complex<From>> & state)
    {
	std::vector<complex<To>> result;
	for (std::size_t i = 0; i < state.size(); i++) {
	    result.push_back(complex<To>(static_cast<To>(state[i].real),
					 static_cast<To>(state[i].imag)));
	}
	return result;
    }

    template<typename To, typename From>
    typename std::enable_if_t<std::is_convertible_v<From, To>,
			      std::vector<To>>
    convertVector(const std::vector<From> & v)
    {
	std::vector<To> result(v.size());
	for (std::size_t i = 0; i < v.size(); i++) {
	    result[i] = static_cast<To>(v[i]);
	}
	return result;
    }


    /**
     * \brief Calculates the hamming weight of an integer.
     * 
     * Counts the number of ones in the binary representation of the integer. 
     */
    unsigned hammingWeight(std::size_t n);

    /**
     * \brief Calculates the binomial coefficient n choose k. 
     */
    unsigned choose(unsigned n, unsigned k);

    /**
     * \brief Generate the lookup vector for a given bit-string
     * length and number of ones.
     *
     * Used in the number preserving simulators
     */
    std::vector<std::size_t> generateLookup(unsigned len, unsigned ones);
    
}


/**
 * \brief Print a std::vector of arbitrary type
 */
template<typename T>
std::ostream& operator<< (std::ostream &stream, const std::vector<T> &vec)
{
    for (unsigned int i = 0; i < vec.size(); i++) {
	stream << vec[i] << std::endl;
    }
    return stream;
}

#endif
