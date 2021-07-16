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
