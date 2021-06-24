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
 * \file default.cpp
 * \brief Contains the implementation of constructors and 
 *        set/get methods of the Default qubits class.
 */


#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <iostream>

template<>
const std::string Qubits<Type::Default, double>::name =
    std::string("Qub<def,double>");

template<>
const std::string Qubits<Type::Default, float>::name =
    std::string("Qub<def,float>");

template<typename Fp>
Qubits<Type::Default, Fp>::Qubits(unsigned nqubits_in) 
    : nqubits{ nqubits_in }, dim{ std::size_t(1) << nqubits },
      state(dim), random(0,1)
{
    // Make the all-zero state
    reset();   
}

template<typename Fp>
Qubits<Type::Default, Fp>::Qubits(const std::vector<complex<Fp>> & state)
    : nqubits{ checkStateSize(state) }, dim{ state.size() }, state{state},
      random(0,1)
{
    //std::cout << "nqubits = " << nqubits << std::endl;
    //std::cout << "dim = " << dim << std::endl;
}

template<typename Fp>
void Qubits<Type::Default, Fp>::reset()
{
    for (std::size_t n=0; n<dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    state[0].real = 1;
}

template<typename Fp>
void Qubits<Type::Default, Fp>::setState(const std::vector<complex<Fp>> & state_in)
{
    if (state_in.size() != dim) {
	std::string msg = "Cannot assign state vector from different ";
	msg += "number of qubits";
	throw std::logic_error(msg);
    }

    // Set the new state vector
    state = state_in;
}

template<typename Fp>
void Qubits<Type::Default, Fp>::setBasisState(std::size_t index)
{
    // Clear the state - set all amplitudes to zero
    for (std::size_t n = 0; n < dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    // Set the amplitude for index to 1
    state[index].real = 1;
}



template<typename Fp>
void Qubits<Type::Default, Fp>::operator = (const Qubits & old)
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
template<typename Fp>
std::vector<complex<Fp>> Qubits<Type::Default, Fp>::getState() const
{
    return state;
}

template<typename Fp>
unsigned Qubits<Type::Default, Fp>::getNumQubits() const
{
    return nqubits;
}


template<typename Fp>
void Qubits<Type::Default, Fp>::print() const
{
    std::cout << "Number of qubits = " << nqubits << std::endl;
    std::cout << state << std::endl;
}


// Explicit instantiations
template class Qubits<Type::Default, float>;
template class Qubits<Type::Default, double>;
