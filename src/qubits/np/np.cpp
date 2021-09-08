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

template<>
const std::string qsl::Qubits<qsl::Type::NP, double>::name =
    std::string("Qub<np,double>");

template<>
const std::string qsl::Qubits<qsl::Type::NP, float>::name =
    std::string("Qub<np,float>");

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::NP, Fp>::Qubits(unsigned nqubits_in, unsigned nones_in) 
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

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::NP, Fp>::Qubits(const std::vector<qsl::complex<Fp>> & state)
    : nqubits{ qsl::checkStateSize(state) }, dim{ state.size() },
      nones{ checkStateNP(state) }, state{ state }, random(0,1)
{
    initLookup();
    //std::cout << "nqubits = " << nqubits << std::endl;
    //std::cout << "dim = " << dim << std::endl;
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::NP, Fp>::reset()
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
void qsl::Qubits<qsl::Type::NP, Fp>::setState(const std::vector<qsl::complex<Fp>> & state_in)
{   
    if (state_in.size() != dim) {
	std::string msg = "Cannot assign state vector from different ";
	msg += "number of qubits";
	throw std::logic_error(msg);
    }

    nones = checkStateNP(state_in);
    initLookup();
    
    // Set the new state vector
    state = state_in;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::NP, Fp>::setBasisState(std::size_t index)
{
    unsigned nones_old = nones;
    nones = qsl::hammingWeight(index);
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
void qsl::Qubits<qsl::Type::NP, Fp>::operator = (const Qubits & old)
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
std::vector<qsl::complex<Fp>> qsl::Qubits<qsl::Type::NP, Fp>::getState() const
{
    return state;
}

template<std::floating_point Fp>
unsigned qsl::Qubits<qsl::Type::NP, Fp>::getNumQubits() const
{
    return nqubits;
}


template<std::floating_point Fp>
unsigned qsl::Qubits<qsl::Type::NP, Fp>::getNumOnes() const
{
    return nones;
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::NP, Fp>::setNumOnes(unsigned nones_in) 
{
    if (nones_in < 1 or nones_in >= nqubits) {
	throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
    }
    
    nones = nones_in;
    initLookup();
    reset();
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::NP, Fp>::initLookup() 
{
    lookup.clear();
    // Generate lookup tables
    // For use in 1-qubit gates and wavefunction collapse
    lookup[{nqubits-1, nones-1}] = qsl::generateLookup(nqubits-1, nones-1);
    // Used for 2-qubit gates
    lookup[{nqubits-2, nones-1}] = qsl::generateLookup(nqubits-2, nones-1);
    lookup[{nqubits-2, nones-2}] = qsl::generateLookup(nqubits-2, nones-2);
    // Used in measurement functions - e.g. generate distribution, collapse
    lookup[{nqubits-1, nones}] = qsl::generateLookup(nqubits-1, nones);
    lookup[{nqubits, nones}] = qsl::generateLookup(nqubits, nones);   
}



template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::NP, Fp>::print(std::ostream & os) const
{
    os << "Number of qubits = " << nqubits << std::endl;
    os << "Number of ones = " << nones << std::endl;
    os << state << std::endl;
}


// Explicit instantiations
template class qsl::Qubits<qsl::Type::NP, float>;
template class qsl::Qubits<qsl::Type::NP, double>;

