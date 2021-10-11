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
 * \file qubits_thread.cpp
 * \brief Contains the implementation of the std::thread-based Qubits class
 *
 */

#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <iostream>
#include <cmath>
#include <string>

template<>
const std::string Qubits<Type::Thread, double>::name =
    std::string("Qub<thread,double>");

template<>
const std::string Qubits<Type::Thread, float>::name =
    std::string("Qub<thread,float>");


template<typename Fp>
void Qubits<Type::Thread, Fp>::pauliX_inner(std::size_t start,
					    std::size_t end,
					    std::size_t k)
{
    for (std::size_t s = 0; s < dim; s += 2*k) { 
	for (std::size_t r = start; r < end; r++) {
	    // Get the indices that need to be switched
	    std::size_t index1 = s + r;
	    std::size_t index2 = s + k + r;
	    std::swap(state[index1], state[index2]);
	}
    }
}

template<typename Fp>
void Qubits<Type::Thread, Fp>::pauliX_0(std::size_t start,
					std::size_t end)
{
    for (std::size_t s = start; s < end; s += 2) { 
	// Get the indices that need to be switched
	std::swap(state[s], state[s+1]);
    }
}
template<typename Fp>
void Qubits<Type::Thread, Fp>::pauliX_1(std::size_t start,
					std::size_t end)
{
    for (std::size_t s = start; s < end; s += 4) { 
	// Get the indices that need to be switched
	std::swap(state[s], state[s+2]);
	std::swap(state[s+1], state[s+3]);
    }
}
template<typename Fp>
void Qubits<Type::Thread, Fp>::pauliX_2(std::size_t start,
					std::size_t end)
{
    for (std::size_t s = start; s < end; s += 8) { 
	// Get the indices that need to be switched
	std::swap(state[s], state[s+4]);
	std::swap(state[s+1], state[s+5]);
	std::swap(state[s+2], state[s+6]);
	std::swap(state[s+3], state[s+7]);
    }
}
template<typename Fp>
void Qubits<Type::Thread, Fp>::pauliX_3(std::size_t start,
					std::size_t end)
{
    for (std::size_t s = start; s < end; s += 16) { 
	// Get the indices that need to be switched
	std::swap(state[s], state[s+8]);
	std::swap(state[s+1], state[s+9]);
	std::swap(state[s+2], state[s+10]);
	std::swap(state[s+3], state[s+11]);
	std::swap(state[s+4], state[s+12]);
	std::swap(state[s+5], state[s+13]);
	std::swap(state[s+6], state[s+14]);
	std::swap(state[s+7], state[s+15]);
    }
}




template<typename Fp>
void Qubits<Type::Thread, Fp>::pauliX_outer(std::size_t start,
					  std::size_t end,
					  std::size_t k)
{
    for (std::size_t s = start; s < end; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the indices that need to be switched
	    std::size_t index1 = s + r;
	    std::size_t index2 = s + k + r;
	    std::swap(state[index1], state[index2]);
	}
    }
}
 
template<typename Fp>
Qubits<Type::Thread, Fp>::Qubits(unsigned nqubits_in) 
    : nqubits{ nqubits_in }, dim{ std::size_t(1) << nqubits },
      state(dim), random(0,1)
{		
    // Make the all-zero state
    reset();   
}

///\todo This is sure to contain bugs! Very hastily written
template<typename Fp>
Qubits<Type::Thread, Fp>::Qubits(const Qubits & old) 
    : nqubits{ old.nqubits }, dim{ std::size_t(1) << nqubits },
      state(old.state), random(0,1)

{ }


template<typename Fp>
Qubits<Type::Thread, Fp>::Qubits(const std::vector<complex<Fp>> & state)
    : nqubits{ checkStateSize(state) }, dim{ state.size() }, state{state},
      random(0,1)
{ }

template<typename Fp>
void Qubits<Type::Thread, Fp>::reset()
{
    for (std::size_t n = 0; n < dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    state[0].real = 1;
}

template<typename Fp>
void Qubits<Type::Thread, Fp>::setState(const std::vector<complex<Fp>> & state_in)
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
void Qubits<Type::Thread, Fp>::setBasisState(std::size_t index)
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
void Qubits<Type::Thread, Fp>::operator = (const Qubits & old) {

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
std::vector<complex<Fp>> Qubits<Type::Thread, Fp>::getState() const
{
    return state;
}

template<typename Fp>
unsigned Qubits<Type::Thread, Fp>::getNumQubits() const
{
    return nqubits;
}


template<typename Fp>
void Qubits<Type::Thread, Fp>::print() const
{
    std::cout << "Number of qubits = " << nqubits << std::endl;
    std::cout << state << std::endl;
}


template<typename Fp>
void Qubits<Type::Thread, Fp>::pauliX(unsigned targ)
{
    // Set up the variables and start the thread going
    std::size_t k = 1 << targ;

    const std::size_t hardware_concurrency = 4;

    if (targ < nqubits-2) {
	std::vector<std::future<void>> results;
	for (std::size_t t = 0; t < hardware_concurrency; t++) {
	    std::size_t start = t*(dim/hardware_concurrency);
	    std::size_t end = (t+1)*(dim/hardware_concurrency);
	    results.push_back(pool.submit(std::bind(&Qubits::pauliX_outer,
	    					    this, start, end, k)));
	}
	for (std::size_t t = 0; t < hardware_concurrency; t++) {
	    results[t].wait();
	}
    } else {

// std::future<void> r1 = pool.submit(std::bind(&Qubits::pauliX_inner,
	// 					     this, 0, k/4, k));
	// std::future<void> r2 = pool.submit(std::bind(&Qubits::pauliX_inner,
	// 					     this, k/4, k/2, k));
	// std::future<void> r3 = pool.submit(std::bind(&Qubits::pauliX_inner,
	// 					     this, k/2, 3*k/4, k));
	// std::future<void> r4 = pool.submit(std::bind(&Qubits::pauliX_inner,
	// 					     this, 3*k/4, k, k));
	// r1.wait();
	// r2.wait();
	// r3.wait();
	// r4.wait();
	pauliX_outer(0,dim,k);
    }
}

template<typename Fp>
void Qubits<Type::Thread, Fp>::controlNot(unsigned ctrl, unsigned targ)
{
    std::size_t small_bit = 1 << std::min(ctrl, targ);
    std::size_t large_bit = 1 << std::max(ctrl, targ);

    std::size_t mid_incr = (small_bit << 1);
    std::size_t high_incr = (large_bit << 1);
    std::size_t targ_bit = (1 << targ);
    std::size_t ctrl_bit = (1 << ctrl);

    // Increment through the indices above largest bit (ctrl or targ)
    for (std::size_t i = 0; i < dim; i += high_incr) {
	// Increment through the middle set of bits
	for (std::size_t j = 0; j < large_bit; j += mid_incr) {
	    // Increment through the low set of bits
	    for (std::size_t k = 0; k < small_bit; k++) {
		// 2x2 matrix multiplication on the zero (i+j+k)
		// and one (i+j+k+targ_bit) indices. 
		//  mat_mul(op, state, i+j+k, i+j+k+targ_bit);
		std::size_t indexUp = i + j + k + ctrl_bit;
		std::size_t indexLo = indexUp + targ_bit;
		
		//complex temp = state[indexUp];
		//state[indexUp] = state[indexLo];
		//state[indexLo] = temp;
		std::swap(state[indexUp], state[indexLo]);
	    }
	}
    }    
}

// Explicit instantiations
template class Qubits<Type::Thread, float>;
template class Qubits<Type::Thread, double>;
