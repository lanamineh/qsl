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
 * \file resize.cpp
 * \brief Contains the implementation of constructors and 
 *        set/get methods of the Resize qubits class.
 */


#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include "qsl/utils/quantum.hpp"
#include <iostream>

template<>
const std::string qsl::Qubits<qsl::Type::Resize, double>::name =
    std::string("Qub<resize,double>");

template<>
const std::string qsl::Qubits<qsl::Type::Resize, float>::name =
    std::string("Qub<resize,float>");

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::Resize, Fp>::Qubits(unsigned nqubits_in) 
    : nqubits{ nqubits_in }, dim{ std::size_t(1) << nqubits },
      state(dim), random(0,1)
{
    // Make the all-zero state
    reset();   
}

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::Resize, Fp>::Qubits(const std::vector<qsl::complex<Fp>> & state)
    : nqubits{ qsl::checkStateSize(state) }, dim{ state.size() }, state{state},
      random(0,1)
{
    //std::cout << "nqubits = " << nqubits << std::endl;
    //std::cout << "dim = " << dim << std::endl;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::reset()
{
    for (std::size_t n=0; n<dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    state[0].real = 1;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::setState(const std::vector<qsl::complex<Fp>> & state_in)
{
    // For the resizeable simulator, you can set a different sized
    // state vector without causing an exception (unless the state is
    // an invalid size)
    dim = state_in.size();
    nqubits = checkStateSize(state_in);
   
    // Set the new state vector
    state = state_in;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::setBasisState(std::size_t index)
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
void qsl::Qubits<qsl::Type::Resize, Fp>::addQubit()
{
    // Push zeros to the back of the state vector
    for (std::size_t k = 0; k < dim; k++) {
	state.push_back({0,0});
    }

    // Store the new number of qubits and the state vector size
    nqubits++;
    dim = state.size();
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::operator = (const Qubits & old)
{
    // For the resizeable simulator, you can assign a different sized
    // simulator object without causing an exception (unless the state is
    // an invalid size)
    dim = old.dim;
    nqubits = old.nqubits;

    // Set the new state vector
    state = old.state;
}

/// Get the state vector associated to the qubits
template<std::floating_point Fp>
std::vector<qsl::complex<Fp>> qsl::Qubits<qsl::Type::Resize, Fp>::getState() const
{
    return state;
}

template<std::floating_point Fp>
unsigned qsl::Qubits<qsl::Type::Resize, Fp>::getNumQubits() const
{
    return nqubits;
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::print(std::ostream & os) const
{
    os << "Number of qubits = " << nqubits << std::endl;
    os << state << std::endl;
}


// Explicit instantiations
template class qsl::Qubits<qsl::Type::Resize, float>;
template class qsl::Qubits<qsl::Type::Resize, double>;
