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
 * \file random.cpp
 * \brief Implementation of random number related functions
 *
 *
 */

#include "qsl/utils/random.hpp"
#include "qsl/utils/quantum.hpp"
#include "qsl/utils/misc.hpp"

namespace qsl {
 
    template<std::floating_point Fp>
    Random<Fp>::Random(Fp a, Fp b)
    {
	///\todo Is std::random_device OK?
	std::random_device r;
	
	generator = std::mt19937(r());
	distribution = std::uniform_real_distribution<Fp>(a,b);
    }

    template<std::floating_point Fp>
    Fp Random<Fp>::getNum() {
	return distribution(generator);
    }


    template<std::floating_point Fp>
    std::vector<complex<Fp>> makeRandomState(unsigned nqubits)
    {
	std::size_t dim = 1 << nqubits;

	qsl::Random<Fp> random(-1,1);

    
    
	std::vector<complex<Fp>> state;
	for(std::size_t i=0; i<dim; i++) {
	    Fp val_real = random.getNum();
	    Fp val_imag = random.getNum();
	    state.push_back(complex(val_real, val_imag));
	}

	// Normalise the state vector
	normalise(state);
    
	return state;
    }

    ///\todo This needs testing
    template<std::floating_point Fp>
    std::vector<complex<Fp>> makeRandomNPState(unsigned nqubits)
    {
	// Make a random number of ones
	std::random_device r;    
	std::mt19937 generator{r()};
	std::uniform_int_distribution<unsigned> distribution(1,nqubits-1);
	unsigned nones = distribution(generator);

	return makeRandomNPState<Fp>(nqubits,nones);
    }

    ///\todo This needs testing
    template<std::floating_point Fp>
    std::vector<complex<Fp>> makeRandomNPState(unsigned nqubits, unsigned nones)
    {
	Random<Fp> random{-1,1};
	
	std::size_t dim = 1 << nqubits;
	std::size_t x = (1ULL << nones) - 1;
	std::size_t end = x << (nqubits - nones);

	// Make state vector with correct length
	std::vector<complex<Fp>> state(dim);
    
	while (x <= end) {

	    Fp val_real = random.getNum();
	    Fp val_imag = random.getNum();
	    state[x] = complex(val_real, val_imag);	
	    next(x);
	}
    
	// Normalise the state vector
	normalise(state);
    
	return state;
    }

    

    template<std::floating_point Fp>
    std::vector<Fp> makeRandomPhases(int test_len)
    {
	std::vector<Fp> phase_list;

	qsl::Random<Fp> random(-M_PI, M_PI);
    
	for (int i = 0; i < test_len; i++) {
	    phase_list.push_back(random.getNum());
	}
	return phase_list;
    }

    /// Explicit instantiations
    template class qsl::Random<float>;
    template class qsl::Random<double>;

    template std::vector<complex<float>> makeRandomState(unsigned nqubits);
    template std::vector<complex<double>> makeRandomState(unsigned nqubits);

    template std::vector<complex<float>> makeRandomNPState(unsigned nqubits);
    template std::vector<complex<double>> makeRandomNPState(unsigned nqubits);

    template std::vector<complex<float>> makeRandomNPState(unsigned nqubits,
							   unsigned nones);
    template std::vector<complex<double>> makeRandomNPState(unsigned nqubits,
							    unsigned nones);

    
    template std::vector<float> makeRandomPhases(int nqubits);
    template std::vector<double> makeRandomPhases(int nqubits);

}
