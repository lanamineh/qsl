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
 * \file np.cpp
 * \brief Contains the implementation of the number preserved Qubits class.
 *
 */


#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include "qsl/utils/quantum.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

namespace qsl {

    template<>
    const std::string Qubits<Type::NP, double>::name =
	std::string("Qub<np,double>");

    template<>
    const std::string Qubits<Type::NP, float>::name =
	std::string("Qub<np,float>");

    template<std::floating_point Fp>
    Qubits<Type::NP, Fp>::Qubits(unsigned nqubits_in, unsigned nones_in) 
	: nqubits{ nqubits_in }, dim{ std::size_t(1) << nqubits },
	  nones{ nones_in }, state(dim), random(0,1)
    {
	/// \todo Figure out better place to put this?
	if (nones < 1 or nones >= nqubits) {
	    throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
	}

	initLookup();
	// Make the computational basis state with the lowest value
	reset();   
    }


    /**
     * \brief Check if the state is number preserving and find the number
     * of ones if it is.
     *
     * This function is probably inefficient - needs improving.
     * Also might make it a member function of the class 
     */
    template<std::floating_point Fp>
    unsigned checkStateNP(const std::vector<complex<Fp>> & state)
    {
	unsigned nones = 0;
	bool found = false;
    
	for (std::size_t i = 0; i < state.size(); i++) {
	    Fp amp = state[i].real * state[i].real +
		state[i].imag * state[i].imag;
	    // If amplitude is non-zero store the number of ones
	    if (amp != 0) {
		unsigned weight = hammingWeight(i);
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



    template<std::floating_point Fp>
    Qubits<Type::NP, Fp>::Qubits(const std::vector<complex<Fp>> & state)
	: nqubits{ checkStateSize(state) }, dim{ state.size() },
	  nones{ checkStateNP(state) }, state{ state }, random(0,1)
    {
	if (nones < 1 or nones >= nqubits) {
	    throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
	}

	initLookup();
	//std::cout << "nqubits = " << nqubits << std::endl;
	//std::cout << "dim = " << dim << std::endl;
    }


    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::reset()
    {
	for (std::size_t n = 0; n < dim; n++) {
	    state[n].real = 0;
	    state[n].imag = 0;
	}
    
	// Starting number (nones 1s in the least significant positions)
	std::size_t idx = (1ULL << nones) - 1ULL;
	state[idx].real = 1;
    }

    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::setState(const std::vector<complex<Fp>> & state_in)
    {   
	if (state_in.size() != dim) {
	    std::string msg = "Cannot assign state vector from different ";
	    msg += "number of qubits";
	    throw std::logic_error(msg);
	}

	nones = checkStateNP(state_in);
	if (nones < 1 or nones >= nqubits) {
	    throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
	}

	initLookup();
    
	// Set the new state vector
	state = state_in;
    }

    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::setBasisState(std::size_t index)
    {
	unsigned nones_old = nones;
	nones = hammingWeight(index);
	if (nones_old != nones) {
	    initLookup();
	}
    
	// Clear the state - set all amplitudes to zero
	for (std::size_t n = 0; n < dim; n++) {
	    state[n].real = 0;
	    state[n].imag = 0;
	}
	// Set the amplitude for index to 1
	state[index].real = 1;
    }



    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::operator = (const Qubits & old)
    {
	// Check if nqubits (and therefore dim) are the same
	if (nqubits != old.nqubits) {
	    std::string msg = "Cannot assign Qubits object vector from different ";
	    msg += "number of qubits";
	    throw std::logic_error(msg);
	}
    
	// Set the new state vector
	state = old.state;
    }

    /// Get the state vector associated to the qubits
    template<std::floating_point Fp>
    std::vector<complex<Fp>> Qubits<Type::NP, Fp>::getState() const
    {
	return state;
    }

    template<std::floating_point Fp>
    unsigned Qubits<Type::NP, Fp>::getNumQubits() const
    {
	return nqubits;
    }


    template<std::floating_point Fp>
    unsigned Qubits<Type::NP, Fp>::getNumOnes() const
    {
	return nones;
    }


    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::setNumOnes(unsigned nones_in) 
    {
	if (nones_in < 1 or nones_in >= nqubits) {
	    throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
	}
    
	nones = nones_in;
	initLookup();
	reset();
    }


    /// Generate the lookup vector for a given bit-string length and number of ones
    [[nodiscard]] std::vector<std::size_t> generateLookup(unsigned len, unsigned ones)
    {
	std::vector<std::size_t> list;
	// Starting number (ones in the least significant positions)
	std::size_t x = (1ULL << ones) - 1ULL;
	// Last number (ones in the most significant positions)
	std::size_t end = x << (len - ones);

	// Loop through all the other numbers
	while(x <= end) {
	    list.push_back(x);
	    next(x);
	}

	return list;
    }


    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::initLookup() 
    {
	lookup.clear();
	// Generate lookup tables
	// For use in 1-qubit gates and wavefunction collapse
	lookup[{nqubits-1, nones-1}] = generateLookup(nqubits-1, nones-1);
	// Used for 2-qubit gates
	lookup[{nqubits-2, nones-1}] = generateLookup(nqubits-2, nones-1);
	lookup[{nqubits-2, nones-2}] = generateLookup(nqubits-2, nones-2);
	// Used in measurement functions - e.g. generate distribution, collapse
	lookup[{nqubits-1, nones}] = generateLookup(nqubits-1, nones);
	lookup[{nqubits, nones}] = generateLookup(nqubits, nones);   
    }



    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::print() const
    {
	std::cout << "Number of qubits = " << nqubits << std::endl;
	std::cout << "Number of ones = " << nones << std::endl;
	std::cout << state << std::endl;
    }


    // Explicit instantiations
    template class Qubits<Type::NP, float>;
    template class Qubits<Type::NP, double>;

}
