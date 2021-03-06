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
 * \brief Contains the implementation of the gates for 
 *        the number preserved Qubits class with OpenMP.
 */


#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <cmath>
#include <algorithm>


/* One-qubit gates ***************************************************/

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::phase(unsigned targ, Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	qsl::complex<Fp> phase{std::cos(angle), std::sin(angle)};

	// Position of target qubit
	std::size_t k = 1 << targ;

	// Create masks for breaking up the number
	std::size_t lower_mask = k - 1;
	std::size_t upper_mask = ~lower_mask;
    
	// Loop through all the other numbers
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-1, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-1, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t upper = x & upper_mask;
	    // Get index for 1 in target position
	    std::size_t index = lower + k + (upper << 1);

	    qsl::complex<Fp> temp = state[index];
	    state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	    state[index].imag = phase.real * temp.imag + phase.imag * temp.real;
	}
    }
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::pauliZ(unsigned targ)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Position of target qubit
	std::size_t k = 1 << targ;

	// Create masks for breaking up the number
	std::size_t lower_mask = k - 1;
	std::size_t upper_mask = ~lower_mask;
    
	// Loop through all the other numbers
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-1, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-1, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t upper = x & upper_mask;
	    // Get index for 1 in target position
	    std::size_t index = lower + k + (upper << 1);

	    state[index].real *= -1;
	    state[index].imag *= -1;
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::rotateZ(unsigned targ, Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	// Position of target qubit
	std::size_t k = 1 << targ;

	// Create masks for breaking up the number
	std::size_t lower_mask = k - 1;
	std::size_t upper_mask = ~lower_mask;

	// Apply e^(-i*angle/2) to |0>
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-1, nones}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-1, nones})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t upper = x & upper_mask;
	    // Get index for 1 in target position
	    std::size_t index = lower + (upper << 1);

	    qsl::complex<Fp> temp = state[index];
	    state[index].real = cos * temp.real + sin * temp.imag;
	    state[index].imag = cos * temp.imag - sin * temp.real;
	}
    }

    // Apply e^(i*angle/2) to |1>
    // Because of the way the NP simulator is coded, this requires
    // another loop, so we may as well use the phase function
    phase(targ, angle/2);
}


/* Two-qubit gates ***************************************************/

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::controlZ(unsigned ctrl, unsigned targ)
{
    // Gate does nothing if there is only one 1.
    if (nones < 2) {
	return;
    }
#pragma omp parallel num_threads(nthreads)
    {
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
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-2}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-2})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit (which represent having 1
	    // on ctrl and targ).
	    std::size_t index = lower + small_bit + (mid << 1) + large_bit + (upper << 2);
	    state[index].real *= -1; 
	    state[index].imag *= -1;
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::controlRotateZ(unsigned ctrl,
						       unsigned targ,
						       Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	// Find the bit positions of ctrl and targ
	std::size_t small_bit = 1 << std::min(ctrl, targ);
	std::size_t large_bit = 1 << std::max(ctrl, targ);
	std::size_t ctrl_bit = 1 << ctrl;
    
	// Create masks for the 3 sections that the bit string will be broken into 
	std::size_t lower_mask = small_bit - 1;
	std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
	std::size_t upper_mask = ~(lower_mask | mid_mask);
        
	// Add the e^{-i*angle/2} phase to |01>
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit (which represent having 1
	    // on ctrl and targ).
	    std::size_t index = lower + (mid << 1) + (upper << 2) + ctrl_bit;
	    qsl::complex<Fp> temp = state[index];
	    state[index].real = cos * temp.real + sin * temp.imag;
	    state[index].imag = cos * temp.imag - sin * temp.real;
	}
    }
    // Now add the e^{i*angle/2} phase to |11>
    // Because of the way the NP simulator is coded, this requires
    // another loop, so we may as well use the controlPhase function
    controlPhase(ctrl, targ, angle/2);
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::controlPhase(unsigned ctrl, unsigned targ,
						     Fp angle)
{
    // Gate does nothing if there is only one 1.
    if (nones < 2) {
	return;
    }

#pragma omp parallel num_threads(nthreads)
    {    
	qsl::complex<Fp> phase{std::cos(angle), std::sin(angle)};

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
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-2}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-2})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit (which represent having 1
	    // on ctrl and targ).
	    std::size_t index = lower + small_bit + (mid << 1) + large_bit + (upper << 2);

	    qsl::complex<Fp> amp = state[index];
	    state[index].real = phase.real * amp.real - phase.imag * amp.imag;
	    state[index].imag = phase.real * amp.imag + phase.imag * amp.real;
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::swap(unsigned q1, unsigned q2)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Find the bit positions
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);

	// Create masks for the 3 sections that the bit string will be broken into 
	std::size_t lower_mask = small_bit - 1;
	std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
	std::size_t upper_mask = ~(lower_mask | mid_mask);
        
	// Loop through all the other numbers with num_ones 1s 
	// then break down that number into 3 parts to go on either side
	// of q1 and q2
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit
	    std::size_t index01 = lower + small_bit + (mid << 1) + (upper << 2);
	    std::size_t index10 = lower + (mid << 1) + large_bit + (upper << 2);

	    std::swap(state[index01], state[index10]);
	}
    }
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::fswap(unsigned q1, unsigned q2)
{
    // Since we would have to use two lookup tables (and therefore two loops)
    // to implement the swapping of |01> and |10> and then the phase shift on |11>,
    // it is best just to use separate functions for them
    swap(q1, q2);
    controlZ(q1, q2);
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::npRotateX(unsigned q1,
						  unsigned q2,
						  Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	// Find the bit positions
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);
	std::size_t q1_bit = 1 << q1;
	std::size_t q2_bit = 1 << q2;
    
	// Create masks for the 3 sections that the bit string will be broken into 
	std::size_t lower_mask = small_bit - 1;
	std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
	std::size_t upper_mask = ~(lower_mask | mid_mask);
        
	// Loop through all the other numbers with num_ones 1s 
	// then break down that number into 3 parts to go on either side
	// of q1 and q2
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit
	    std::size_t index1 = lower + (mid << 1) + (upper << 2) + q1_bit;
	    std::size_t index2 = lower + (mid << 1) + (upper << 2) + q2_bit;

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

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::npRotateY(unsigned q1,
						  unsigned q2,
						  Fp angle)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Store variables
	Fp cos = std::cos(angle/2);
	Fp sin = std::sin(angle/2);
    
	// Find the bit positions
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);
	std::size_t q1_bit = 1 << q1;
	std::size_t q2_bit = 1 << q2;
    
	// Create masks for the 3 sections that the bit string will be broken into 
	std::size_t lower_mask = small_bit - 1;
	std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
	std::size_t upper_mask = ~(lower_mask | mid_mask);
        
	// Loop through all the other numbers with num_ones 1s 
	// then break down that number into 3 parts to go on either side
	// of q1 and q2
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit
	    std::size_t index1 = lower + (mid << 1) + (upper << 2) + q1_bit;
	    std::size_t index2 = lower + (mid << 1) + (upper << 2) + q2_bit;

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

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::npHadamard(unsigned q1,
						   unsigned q2)
{
#pragma omp parallel num_threads(nthreads)
    {

	// Find the bit positions
	std::size_t small_bit = 1 << std::min(q1, q2);
	std::size_t large_bit = 1 << std::max(q1, q2);
	std::size_t q1_bit = 1 << q1;
	std::size_t q2_bit = 1 << q2;
    
	// Create masks for the 3 sections that the bit string will be broken into 
	std::size_t lower_mask = small_bit - 1;
	std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
	std::size_t upper_mask = ~(lower_mask | mid_mask);
        
	// Loop through all the other numbers with num_ones 1s 
	// then break down that number into 3 parts to go on either side
	// of q1 and q2
#pragma omp for
	for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-1}).size(); i++) {
	    std::size_t x = lookup.at({nqubits-2, nones-1})[i];
	    std::size_t lower = x & lower_mask;
	    std::size_t mid = x & mid_mask;
	    std::size_t upper = x & upper_mask;

	    // Calculate the index by adding together the 3 shifted sections
	    // of x and small_bit and large_bit
	    std::size_t index1 = lower + (mid << 1) + (upper << 2) + q1_bit;
	    std::size_t index2 = lower + (mid << 1) + (upper << 2) + q2_bit;

	    // Store the values of the amplitudes
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


// Explicit instantiations
template class qsl::Qubits<qsl::Type::OmpNP, float>;
template class qsl::Qubits<qsl::Type::OmpNP, double>;

