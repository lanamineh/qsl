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
 * \brief Contains the implementation of gates for 
 *        the OpenMP-based Qubits class.
 */

#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <cmath>
#include <omp.h>
#include <mutex>


std::mutex mtx;           // mutex for critical section

/* One-qubit gates ***************************************************/

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::pauliX(unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {
	std::size_t k = 1 << targ;

	if(targ < nqubits-2) {	
#pragma omp for
	    for (std::size_t s = 0; s < dim; s += 2*k) { 
		for (std::size_t r = 0; r < k; r++) {
		    // Get the indices that need to be switched
		    std::size_t index1 = s + r;
		    std::size_t index2 = s + k + r;
	    
		    std::swap(state[index1], state[index2]);
		}
	    }
	} else {
	    for (std::size_t s = 0; s < dim; s += 2*k) { 
#pragma omp for
		for (std::size_t r = 0; r < k; r++) {
		    // Get the indices that need to be switched
		    std::size_t index1 = s + r;
		    std::size_t index2 = s + k + r;
		    
		    std::swap(state[index1], state[index2]);
		}
	    }
	}
    }
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::phase(unsigned targ, Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {
	qsl::complex<Fp> phase{std::cos(angle), std::sin(angle)};   
	std::size_t k = 1 << targ;
#pragma omp for
	for (std::size_t s = 0; s < dim; s += 2*k) { 
	    for (std::size_t r = 0; r < k; r++) {
		// Get the index of |1>
		std::size_t index = s + k + r;
	    
		//state[index] *= phase;
		qsl::complex<Fp> temp = state[index];
		state[index].real
		    = phase.real * temp.real - phase.imag * temp.imag;
		state[index].imag
		    = phase.real * temp.imag + phase.imag * temp.real;
	    }
	}
    }
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::rotateX(unsigned targ, Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {
	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
	
	std::size_t k = 1 << targ;
	// Replaced std::size_t with std::size_t for MSCV openmp
	///\todo Check the omp standard for for loop spec
#pragma omp for
	for (std::size_t s = 0; s < dim; s += 2*k) { 
	    for (std::size_t r = 0; r < k; r++) {

		// Get the index of |0> and |1>
		std::size_t index_0 = s + r;
		std::size_t index_1 = s + k + r;

		// Store the values of |0> and |1> amplitudes
		qsl::complex<Fp> a0 = state[index_0];
		qsl::complex<Fp> a1 = state[index_1];

		// Write the new |0> amplitude
		state[index_0].real = a0.real * cos + a1.imag * sin;
		state[index_0].imag = a0.imag * cos - a1.real * sin;

		// Write the new |1> amplitude
		state[index_1].real = a1.real * cos + a0.imag * sin;
		state[index_1].imag = a1.imag * cos - a0.real * sin;
	    
	    }
	}
    }
}

/* Two-qubit gates ***************************************************/

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::controlNot(unsigned ctrl, unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {	
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t targ_bit = (1 << targ);
	std::size_t ctrl_bit = (1 << ctrl);

#pragma omp for
	// Increment through the indices above largest bit (ctrl or targ)
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |11> indices
		    std::size_t index1 = i + j + k + ctrl_bit;
		    std::size_t index2 = index1 + targ_bit;

		    std::swap(state[index1], state[index2]);
		}
	    }
	}
    }    
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::controlPhase(unsigned ctrl,
						   unsigned targ,
						   Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {
	qsl::complex phase{std::cos(angle), std::sin(angle)};
    
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);

	std::size_t outcome = (1 << targ) + (1 << ctrl);
    
#pragma omp for
	// Increment through the indices above largest bit (ctrl or targ)
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    std::size_t index = i + j + k + outcome;
		    // state[index] *= phase;
		    qsl::complex<Fp> amp = state[index];
		    state[index].real
			= phase.real * amp.real - phase.imag * amp.imag;
		    state[index].imag
			= phase.real * amp.imag + phase.imag * amp.real;
		}
	    }
	}    
    }
}


// Explicit instantiations
template class qsl::Qubits<qsl::Type::Omp, float>;
template class qsl::Qubits<qsl::Type::Omp, double>;

