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
 * \file quantum.cpp
 * \brief Implementation of functions in quantum.hpp
 *
 */

#include <random>
#include <vector>
#include <string>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/misc.hpp"

namespace qsl {

    template<std::floating_point Fp = double>
    complex<Fp> innerProduct(const std::vector<complex<Fp>> & v,
			     const std::vector<complex<Fp>> & w)
    {
	if(v.size() != w.size()) {
	    std::string msg = "Cannot compute inner product between vectors ";
	    msg += "of different sizes: ";
	    msg += std::to_string(v.size()) + " and " + std::to_string(w.size());
	    throw std::logic_error(msg);
	}

	complex<Fp> result{ 0,0 };
	for(std::size_t n=0; n<v.size(); n++) {
	    result.real += v[n].real * w[n].real + v[n].imag * w[n].imag;
	    result.imag += v[n].real * w[n].imag - v[n].imag * w[n].real;
	}

	return result;
    }

    template<std::floating_point Fp = double>
    Fp norm(const std::vector<complex<Fp>> & v)
    {
	// Find the norm of the vector
	Fp inner_prod = 0;
	for (std::size_t i = 0; i < v.size(); i++) {
	    inner_prod += v[i].real * v[i].real + v[i].imag * v[i].imag;
	}

	return std::sqrt(inner_prod);
    
    }

    /**
     * \brief Compute arccos after checking argument
     *
     * The function is necessary in case a number rounds 
     * to something bigger than 1 or smaller than -1
     *
     */
    double acosSafe(double value)
    {
	if (value <= -1.0) {
	    return M_PI;
	} else if (value >= 1.0) {
	    return 0;
	} else {
	    return std::acos(value);
	}
    }

    template<std::floating_point Fp = double>
    Fp fubiniStudy(const std::vector<complex<Fp>> & v,
		   const std::vector<complex<Fp>> & w)
    {
	Fp numerator = abs(innerProduct(v,w));
	Fp denominator = norm(v) * norm(w); 
	return acosSafe(numerator/denominator);
    }

    template<std::floating_point Fp = double>
    Fp normalise(std::vector<complex<Fp>> &state)
    {
	// Find the norm of the vector
	Fp factor = norm(state);
   
	// Divide by the norm;
	for (std::size_t i = 0; i < state.size(); i++) {
	    state[i].real /= factor;
	    state[i].imag /= factor;
	}

	return factor;
    }

    template<std::floating_point Fp>
    unsigned checkStateSize(const std::vector<complex<Fp>> & state) {

	// Counts the number of ones in the binary representation of dim
	// If this number is not 1, then dim is not a power of 2 and
	// is therefore invalid
	unsigned one_count = 0;
	// Counts the length of the binary representation of dim
	unsigned len = 0;
    
	std::size_t n = state.size();
	while(n > 0) {
	    len += 1;
	    one_count += n & 1;
	    n >>= 1;
	}
    
	if(one_count != 1) {
	    throw std::logic_error("State vector size is not a power of two");
	}

	return len - 1;

    }

    template<std::floating_point Fp>
    unsigned checkStateNP(const std::vector<qsl::complex<Fp>> & state)
    {
	unsigned nones = 0;
	bool found = false;
    
	for (std::size_t i = 0; i < state.size(); i++) {
	    Fp amp = state[i].real * state[i].real +
		state[i].imag * state[i].imag;
	    // If amplitude is non-zero store the number of ones
	    if (amp != 0) {
		unsigned weight = qsl::hammingWeight(i);
		if (found == false) {
		    nones = weight;
		    found = true;
		}
		else if (nones != weight) {
		    throw std::logic_error("Input state is not number preserving.");
		}
	    }
	}

	return nones;
    }

    /// Explicit instantiations
    template float fubiniStudy(const std::vector<complex<float>> & v,
			       const std::vector<complex<float>> & w);
    template double fubiniStudy(const std::vector<complex<double>> & v,
				const std::vector<complex<double>> & w);

    template float normalise(std::vector<complex<float>> &state);
    template double normalise(std::vector<complex<double>> &state);

    template unsigned checkStateSize(const std::vector<complex<float>> & state);
    template unsigned checkStateSize(const std::vector<complex<double>> & state);

    template unsigned checkStateNP(const std::vector<complex<float>> & state);
    template unsigned checkStateNP(const std::vector<complex<double>> & state);

    template float norm(const std::vector<complex<float>> & v);
    template double norm(const std::vector<complex<double>> & v);

    
}
