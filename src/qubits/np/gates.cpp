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
 * \file gates.cpp
 * \brief Contains the implementation of the gates for 
 *        the number preserved Qubits class.
 */


#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <cmath>
#include <algorithm>

namespace qsl {

    /* One-qubit gates ***************************************************/

    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::phase(unsigned targ, Fp angle)
    {
	complex<Fp> phase = complex(std::cos(angle), std::sin(angle));

	// Position of target qubit
	std::size_t k = 1 << targ;

	// Create masks for breaking up the number
	std::size_t lower_mask = k - 1;
	std::size_t upper_mask = ~lower_mask;
    
	// Loop through all the other numbers
	//#pragma omp parallel for
	for (std::size_t i = 0; i < lookup.at({nqubits-1, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-1, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t upper = x & upper_mask;
	    // Get index for 1 in target position
	    std::size_t index = lower + k + (upper << 1);

	    complex temp = state[index];
	    state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	    state[index].imag = phase.real * temp.imag + phase.imag * temp.real;
	}
    }

    /* Two-qubit gates ***************************************************/

    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::controlPhase(unsigned ctrl, unsigned targ,
					    Fp angle)
    {
	// Gate does nothing if there is only one 1.
	if (nones < 2) {
	    return;
	}
    
	complex<Fp> phase = complex(std::cos(angle), std::sin(angle));

	// Find the bit positions of ctrl and targ
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	// Create masks for the 3 sections that the bit string will be broken into 
	std::size_t lower_mask = small_bit - 1;
	std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
	std::size_t upper_mask = ~(lower_mask | mid_mask);
        
	// Loop through all the other numbers with num_ones 1s 
	// then break down that number into 3 parts to go on either side
	// of q1 and q2
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-2}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-2})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit (which represent having 1
	    // on ctrl and targ).
	    std::size_t index = lower + small_bit + (mid << 1) + large_bit + (upper << 2);

	    complex amp = state[index];
	    state[index].real = phase.real * amp.real - phase.imag * amp.imag;
	    state[index].imag = phase.real * amp.imag + phase.imag * amp.real;
	}
    }


    template<std::floating_point Fp>
    void Qubits<Type::NP, Fp>::swap(unsigned q1, unsigned q2)
    {
	// Find the bit positions of ctrl and targ
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);

	// Create masks for the 3 sections that the bit string will be broken into 
	std::size_t lower_mask = small_bit - 1;
	std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
	std::size_t upper_mask = ~(lower_mask | mid_mask);
        
	// Loop through all the other numbers with num_ones 1s 
	// then break down that number into 3 parts to go on either side
	// of q1 and q2
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit (which represent having 1
	    // on ctrl and targ).
	    std::size_t index01 = lower + small_bit + (mid << 1) + (upper << 2);
	    std::size_t index10 = lower + (mid << 1) + large_bit + (upper << 2);

	    std::swap(state[index01], state[index10]);
	}
    }


    // Explicit instantiations
    template class Qubits<Type::NP, float>;
    template class Qubits<Type::NP, double>;

}
