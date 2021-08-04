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
 * \file stategen.hpp
 * \brief State generators for use with the Verify class
 */

#ifndef QSL_STATEGEN_HPP
#define QSL_STATEGEN_HPP

namespace qsl {

    /**
     * \brief State generator object object concept
     *
     * All state generators should conform to this specification.
     * A state generator must
     *
     *   1) contain a default constructor which initialises an internal 
     *      state
     *   2) have a getState (const) method which returns the current internal
     *      state.
     * 
     * The state generator can optionally have a configureState method 
     * (with any parameter list),  which generators a new internal state
     *
     * \todo Need to add something about the floating point type it 
     * returns, so as to make it compatible with the simulators
     */
    template<typename T>
    concept StateGenerator = std::is_default_constructible<T>::value
    	&& requires(T t)
    {
    	// Required member functions
    	t.getState();
    };

    
    /**
     * \brief State generator for number preserved states
     *
     * This class is intended to be used as the InitState parameter in the
     * Verify class. It generates a random number preserved state with a 
     * given number of qubits and ones. This state is then used to as the input
     * to the checking algorithms (to initialise the simulators).
     *
     */
    template<std::floating_point Fp>
    class NPStateGen
    {
	std::vector<complex<Fp>> state;
	qsl::Random<Fp> random;
    
    public:

	NPStateGen() : random(-1,1) {}

	void configureState(unsigned nqubits, unsigned nones) {
	    std::cout << "Setting number of ones = " << nones << std::endl;

	    std::size_t dim = 1 << nqubits;

	    std::size_t x = (1ULL << nones) - 1;
	    std::size_t end = x << (nqubits - nones);

	    // Make state vector with correct length
	    state.clear();
	    state.resize(dim);

	    std::cout << "Generating a random state vector" << std::endl;
	    while (x <= end) {

		Fp val_real = random.getNum();
		Fp val_imag = random.getNum();
		state[x] = complex(val_real, val_imag);	
		next(x);
	    }
    
	    // Normalise the state vector
	    normalise(state);
	}

	[[nodiscard]] std::vector<complex<Fp>> getState() const {
	    return state;
	}

    };

    template<std::floating_point Fp>
    class DefaultStateGen
    {
	std::vector<complex<Fp>> state;
	qsl::Random<Fp> random;
    
    public:

	DefaultStateGen() : random(-1,1) {}

	void configureState(unsigned nqubits) {

	    std::size_t dim = 1 << nqubits;

	    // Make state vector with correct length
	    state.clear();

	    std::cout << "Generating a random state vector" << std::endl;
	    for (std::size_t n = 0; n < dim; n++) {
		Fp val_real = random.getNum();
		Fp val_imag = random.getNum();
		state.push_back(complex(val_real, val_imag));	
	    }
    
	    // Normalise the state vector
	    normalise(state);
	}

	[[nodiscard]] std::vector<complex<Fp>> getState() const {
	    return state;
	}

    };

}

#endif
