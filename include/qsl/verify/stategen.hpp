/* 
 * Copyright (C) 2020 Lana Mineh and John Scott.
 *
 * This file is part of QSL, the quantum computer simulator.
 *
 * QSL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QSL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QSL.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * \brief State generators for use with the Verify class
 *
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
