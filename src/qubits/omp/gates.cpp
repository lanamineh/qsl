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
void qsl::Qubits<qsl::Type::Omp, Fp>::pauliY(unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {    
	std::size_t k = 1 << targ;
#pragma omp for
	for (std::size_t s = 0; s < dim; s += 2*k) { 
	    for (std::size_t r = 0; r < k; r++) {
		// Get the indices that need to be modified
		std::size_t index1 = s + r;
		std::size_t index2 = s + k + r;

		qsl::complex<Fp> temp1 = state[index1];
		qsl::complex<Fp> temp2 = state[index2];

		state[index1].real = temp2.imag;
		state[index1].imag = -temp2.real;
		state[index2].real = -temp1.imag;
		state[index2].imag = temp1.real;
	    }
	}
    } 
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::pauliZ(unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {
	std::size_t k = 1 << targ;
#pragma omp for
	for (std::size_t s = 0; s < dim; s += 2*k) { 
	    for (std::size_t r = 0; r < k; r++) {
		// Get the index of |1>
		std::size_t index = s + k + r;
		// Apply -1 phase
		state[index].real *= -1;
		state[index].imag *= -1;
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::hadamard(unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {
	std::size_t k = 1 << targ;
#pragma omp for
	for (std::size_t s = 0; s < dim; s += 2*k) { 
	    for (std::size_t r = 0; r < k; r++) {
		// Get the indices that need to be modified
		const std::size_t index1 = s + r;
		const std::size_t index2 = s + k + r;
		const qsl::complex<Fp> temp1 = state[index1];
		const qsl::complex<Fp> temp2 = state[index2];
		constexpr Fp sqrt2 = std::sqrt(2); 
		state[index1].real = (temp1.real + temp2.real)/sqrt2;
		state[index1].imag = (temp1.imag + temp2.imag)/sqrt2;
		state[index2].real = (temp1.real - temp2.real)/sqrt2;
		state[index2].imag = (temp1.imag - temp2.imag)/sqrt2;
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

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::rotateY(unsigned targ, Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {    
	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
	std::size_t k = 1 << targ;
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
		state[index_0].real = a0.real * cos - a1.real * sin;
		state[index_0].imag = a0.imag * cos - a1.imag * sin;

		// Write the new |1> amplitude
		state[index_1].real = a0.real * sin + a1.real * cos;
		state[index_1].imag = a0.imag * sin + a1.imag * cos;
	    
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::rotateZ(unsigned targ, Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {    
	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
	std::size_t k = 1 << targ;
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
		state[index_0].real = a0.real * cos + a0.imag * sin;
		state[index_0].imag = a0.imag * cos - a0.real * sin;

		// Write the new |1> amplitude
		state[index_1].real = a1.real * cos - a1.imag * sin;
		state[index_1].imag = a1.imag * cos + a1.real * sin;
	    
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
void qsl::Qubits<qsl::Type::Omp, Fp>::controlY(unsigned ctrl,
					       unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {	
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t targ_bit = (1 << targ);
	std::size_t ctrl_bit = (1 << ctrl);
    
	// Increment through the indices above largest bit (ctrl or targ)
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |11> indices
		    std::size_t index1 = i + j + k + ctrl_bit;
		    std::size_t index2 = index1 + targ_bit;

		    qsl::complex<Fp> temp1 = state[index1];
		    qsl::complex<Fp> temp2 = state[index2];
	    
		    state[index1].real = temp2.imag;
		    state[index1].imag = -temp2.real;
		    state[index2].real = -temp1.imag;
		    state[index2].imag = temp1.real;
		}
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::controlZ(unsigned ctrl,
					       unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {
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
		    state[index].real *= -1;
		    state[index].imag *= -1;
		}
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::controlRotateX(unsigned ctrl,
						     unsigned targ,
						     Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t targ_bit = (1 << targ);
	std::size_t ctrl_bit = (1 << ctrl);

	// Increment through the indices above largest bit (ctrl or targ)
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |11> indices
		    std::size_t index1 = i + j + k + ctrl_bit;
		    std::size_t index2 = index1 + targ_bit;

		    // Store the values of the amplitudes
		    qsl::complex<Fp> a0 = state[index1];
		    qsl::complex<Fp> a1 = state[index2];

		    // Write the new |01> amplitude
		    state[index1].real = a0.real * cos + a1.imag * sin;
		    state[index1].imag = a0.imag * cos - a1.real * sin;

		    // Write the new |11> amplitude
		    state[index2].real = a1.real * cos + a0.imag * sin;
		    state[index2].imag = a1.imag * cos - a0.real * sin;
		}
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::controlRotateY(unsigned ctrl,
						     unsigned targ,
						     Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t targ_bit = (1 << targ);
	std::size_t ctrl_bit = (1 << ctrl);

	// Increment through the indices above largest bit (ctrl or targ)
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |11> indices
		    std::size_t index1 = i + j + k + ctrl_bit;
		    std::size_t index2 = index1 + targ_bit;

		    // Store the values of the amplitudes
		    qsl::complex<Fp> a0 = state[index1];
		    qsl::complex<Fp> a1 = state[index2];

		    // Write the new |01> amplitude
		    state[index1].real = a0.real * cos - a1.real * sin;
		    state[index1].imag = a0.imag * cos - a1.imag * sin;

		    // Write the new |11> amplitude
		    state[index2].real = a0.real * sin + a1.real * cos;
		    state[index2].imag = a0.imag * sin + a1.imag * cos;
		}
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::controlRotateZ(unsigned ctrl,
						     unsigned targ,
						     Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t targ_bit = (1 << targ);
	std::size_t ctrl_bit = (1 << ctrl);

	// Increment through the indices above largest bit (ctrl or targ)
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |11> indices
		    std::size_t index1 = i + j + k + ctrl_bit;
		    std::size_t index2 = index1 + targ_bit;

		    // Store the values of the amplitudes
		    qsl::complex<Fp> a0 = state[index1];
		    qsl::complex<Fp> a1 = state[index2];

		    // Write the new |01> amplitude
		    state[index1].real = a0.real * cos + a0.imag * sin;
		    state[index1].imag = a0.imag * cos - a0.real * sin;

		    // Write the new |11> amplitude
		    state[index2].real = a1.real * cos - a1.imag * sin;
		    state[index2].imag = a1.imag * cos + a1.real * sin;
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

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::controlHadamard(unsigned ctrl,
						      unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t targ_bit = (1 << targ);
	std::size_t ctrl_bit = (1 << ctrl);

	// Increment through the indices above largest bit (ctrl or targ)
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |11> indices
		    std::size_t index1 = i + j + k + ctrl_bit;
		    std::size_t index2 = index1 + targ_bit;

		    const qsl::complex<Fp> temp1 = state[index1];
		    const qsl::complex<Fp> temp2 = state[index2];
		    constexpr Fp sqrt2 = std::sqrt(2); 
		    state[index1].real = (temp1.real + temp2.real)/sqrt2;
		    state[index1].imag = (temp1.imag + temp2.imag)/sqrt2;
		    state[index2].real = (temp1.real - temp2.real)/sqrt2;
		    state[index2].imag = (temp1.imag - temp2.imag)/sqrt2;
		}
	    }
	}
    }
}



template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::swap(unsigned q1, unsigned q2)
{
#pragma omp parallel num_threads(nthreads)
    {

	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t q1_bit = (1 << q1);
	std::size_t q2_bit = (1 << q2);

#pragma omp for
	// Increment through the indices above largest bit 
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |10> indices
		    std::size_t index1 = i + j + k + q1_bit;
		    std::size_t index2 = i + j + k + q2_bit;

		    std::swap(state[index1], state[index2]);
		}
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::fswap(unsigned q1, unsigned q2)
{
#pragma omp parallel num_threads(nthreads)
    {
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t q1_bit = (1 << q1);
	std::size_t q2_bit = (1 << q2);

	// Increment through the indices above largest bit
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |10> indices and swap them
		    std::size_t index1 = i + j + k + q1_bit;
		    std::size_t index2 = i + j + k + q2_bit;

		    std::swap(state[index1], state[index2]);

		    // Get the |11> index and apply a -1 phase
		    std::size_t index = index1 + q2_bit;
		    state[index].real *= -1;
		    state[index].imag *= -1;
		}
	    }
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::npRotateX(unsigned q1,
						unsigned q2,
						Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t q1_bit = (1 << q1);
	std::size_t q2_bit = (1 << q2);

	// Increment through the indices above largest bit 
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |10> indices
		    std::size_t index1 = i + j + k + q1_bit;
		    std::size_t index2 = i + j + k + q2_bit;

		    // Store the values of the amplitudes
		    qsl::complex<Fp> a0 = state[index1];
		    qsl::complex<Fp> a1 = state[index2];

		    // Write the new |01> amplitude
		    state[index1].real = a0.real * cos + a1.imag * sin;
		    state[index1].imag = a0.imag * cos - a1.real * sin;

		    // Write the new |10> amplitude
		    state[index2].real = a1.real * cos + a0.imag * sin;
		    state[index2].imag = a1.imag * cos - a0.real * sin;
		}
	    }
	}
    }
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::npRotateY(unsigned q1,
						unsigned q2,
						Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t q1_bit = (1 << q1);
	std::size_t q2_bit = (1 << q2);

	// Increment through the indices above largest bit 
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |10> indices
		    std::size_t index1 = i + j + k + q1_bit;
		    std::size_t index2 = i + j + k + q2_bit;

		    // Store the values of the amplitudes
		    qsl::complex<Fp> a0 = state[index1];
		    qsl::complex<Fp> a1 = state[index2];

		    // Write the new |01> amplitude
		    state[index1].real = a0.real * cos - a1.real * sin;
		    state[index1].imag = a0.imag * cos - a1.imag * sin;

		    // Write the new |10> amplitude
		    state[index2].real = a0.real * sin + a1.real * cos;
		    state[index2].imag = a0.imag * sin + a1.imag * cos;
		}
	    }
	}
    }
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::npHadamard(unsigned q1,
						 unsigned q2)
{
#pragma omp parallel num_threads(nthreads)
    {

	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);

	std::size_t mid_incr = (small_bit << 1);
	std::size_t high_incr = (large_bit << 1);
	std::size_t q1_bit = (1 << q1);
	std::size_t q2_bit = (1 << q2);

	// Increment through the indices above largest bit 
#pragma omp for
	for (std::size_t i = 0; i < dim; i += high_incr) {
	    // Increment through the middle set of bits
	    for (std::size_t j = 0; j < large_bit; j += mid_incr) {
		// Increment through the low set of bits
		for (std::size_t k = 0; k < small_bit; k++) {
		    // Get the |01> and |10> indices
		    std::size_t index1 = i + j + k + q1_bit;
		    std::size_t index2 = i + j + k + q2_bit;

		    // Store the values of the amplitudes
		    qsl::complex<Fp> a0 = state[index1];
		    qsl::complex<Fp> a1 = state[index2];

		    const qsl::complex<Fp> temp1 = state[index1];
		    const qsl::complex<Fp> temp2 = state[index2];
		    constexpr Fp sqrt2 = std::sqrt(2); 
		    state[index1].real = (temp1.real + temp2.real)/sqrt2;
		    state[index1].imag = (temp1.imag + temp2.imag)/sqrt2;
		    state[index2].real = (temp1.real - temp2.real)/sqrt2;
		    state[index2].imag = (temp1.imag - temp2.imag)/sqrt2;
		}
	    }
	}
    }
}


// Explicit instantiations
template class qsl::Qubits<qsl::Type::Omp, float>;
template class qsl::Qubits<qsl::Type::Omp, double>;

