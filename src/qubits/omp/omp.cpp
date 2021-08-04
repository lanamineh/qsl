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
 * \file omp.cpp
 * \brief Contains the implementation of the OpenMP-based Qubits class
 *
 */

#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include "qsl/utils/quantum.hpp"
#include <iostream>
#include <string>
#include <omp.h>

template<>
const std::string qsl::Qubits<qsl::Type::Omp, double>::name =
    std::string("Qub<omp,double>");

template<>
const std::string qsl::Qubits<qsl::Type::Omp, float>::name =
    std::string("Qub<omp,float>");

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::Omp, Fp>::Qubits(unsigned nqubits_in, unsigned nthreads_in) 
    : nqubits{ nqubits_in }, dim{ std::size_t(1) << nqubits },
      nthreads{ nthreads_in }, state(dim), random(0,1)
{
    // Make the all-zero state
    reset();
}

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::Omp, Fp>::Qubits(const std::vector<qsl::complex<Fp>> & state,
					unsigned nthreads_in)
    : nqubits{ qsl::checkStateSize(state) }, dim{ state.size() },
      nthreads{ nthreads_in }, state{state}, random(0,1)
{
    //std::cout << "nqubits = " << nqubits << std::endl;
    //std::cout << "dim = " << dim << std::endl;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::reset()
{
    for (std::size_t n = 0; n < dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    state[0].real = 1;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::setState(const std::vector<qsl::complex<Fp>> & state_in)
{
    if (state_in.size() != dim) {
	std::string msg = "Cannot assign state vector from different ";
	msg += "number of qubits";
	throw std::logic_error(msg);
    }

    // Set the new state vector
    state = state_in;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::setBasisState(std::size_t index)
{
    // Clear the state - set all amplitudes to zero
    for (std::size_t n = 0; n < dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    // Set the amplitude for index to 1
    state[index].real = 1;
}



template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::operator = (const Qubits & old) {

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
std::vector<qsl::complex<Fp>> qsl::Qubits<qsl::Type::Omp, Fp>::getState() const
{
    return state;
}

template<std::floating_point Fp>
unsigned qsl::Qubits<qsl::Type::Omp, Fp>::getNumQubits() const
{
    return nqubits;
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Omp, Fp>::print(std::ostream & os) const
{
    os << "Number of qubits = " << nqubits << std::endl;
    os << state << std::endl;
}

// Explicit instantiations
template class qsl::Qubits<qsl::Type::Omp, float>;
template class qsl::Qubits<qsl::Type::Omp, double>;

