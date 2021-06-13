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
 * \file qubits_default.cpp
 * \brief Contains the implementation of the Qubits class
 *
 */

#include <ciso646> // For MSVC, to use 'and'/'or' keywords

#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

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


template<typename Fp>
void Qubits<Type::Default, Fp>::pauliX(unsigned targ)
{
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < dim; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the indices that need to be switched
	    std::size_t index1 = s + r;
	    std::size_t index2 = s + k + r;
	    std::swap(state[index1], state[index2]);
	}
    }
}

template<typename Fp>
void Qubits<Type::Default, Fp>::phase(unsigned targ, Fp angle)
{
    complex<Fp> phase = complex(std::cos(angle), std::sin(angle));
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < dim; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the index of |1>
	    std::size_t index = s + k + r;
	    //state[index] *= phase;
	    complex temp = state[index];
	    state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	    state[index].imag = phase.real * temp.imag + phase.imag * temp.real;
	}
    }
}

template<typename Fp>
void Qubits<Type::Default, Fp>::rotateX(unsigned targ, Fp angle)
{
    // Store variables
    Fp cos = std::cos(angle/2);
    Fp sin = std::sin(angle/2);
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < dim; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {

	    // Get the index of |0> and |1>
	    std::size_t index_0 = s + r;
	    std::size_t index_1 = s + k + r;

	    // Store the values of |0> and |1> amplitudes
	    complex<Fp> a0 = state[index_0];
	    complex<Fp> a1 = state[index_1];

	    // Write the new |0> amplitude
	    state[index_0].real = a0.real * cos + a1.imag * sin;
	    state[index_0].imag = a0.imag * cos - a1.real * sin;

	    // Write the new |1> amplitude
	    state[index_1].real = a1.real * cos + a0.imag * sin;
	    state[index_1].imag = a1.imag * cos - a0.real * sin;
	    
	}
    }
}

template<typename Fp>
void Qubits<Type::Default, Fp>::controlNot(unsigned ctrl, unsigned targ)
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

template<typename Fp>
void Qubits<Type::Default, Fp>::controlPhase(unsigned ctrl,
					     unsigned targ,
					     Fp angle)
{
    complex phase = complex(std::cos(angle), std::sin(angle));
    
    std::size_t small_bit = 1 << std::min(ctrl, targ);
    std::size_t large_bit = 1 << std::max(ctrl, targ);

    std::size_t mid_incr = (small_bit << 1);
    std::size_t high_incr = (large_bit << 1);

    std::size_t outcome = (1 << targ) + (1 << ctrl);
    
    // Increment through the indices above largest bit (ctrl or targ)
    for (std::size_t i = 0; i < dim; i += high_incr) {
	// Increment through the middle set of bits
	for (std::size_t j = 0; j < large_bit; j += mid_incr) {
	    // Increment through the low set of bits
	    for (std::size_t k = 0; k < small_bit; k++) {
		// 2x2 matrix multiplication on the zero (i+j+k)
		// and one (i+j+k+targ_bit) indices. 
		//  mat_mul(op, state, i+j+k, i+j+k+targ_bit);
		std::size_t index = i + j + k + outcome;

		complex amp = state[index];
		state[index].real = phase.real * amp.real - phase.imag * amp.imag;
		state[index].imag = phase.real * amp.imag + phase.imag * amp.real;
	    }
	}
    }    
}


template<typename Fp>
void Qubits<Type::Default, Fp>::collapse(unsigned targ, unsigned outcome,
					 Fp factor)
{
    // Collapse to outcome and renormalise the state vector simultaneously   
    std::size_t k = 1 << targ;
    std::size_t res = outcome * k;
    std::size_t not_res = (outcome ^ 1) * k; 

    // Loop through the state vector and renormalise amplitudes associated with
    // the correct outcome and zero out those with the opposite outcome.
    for (std::size_t s = 0; s < dim; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the indices that need to be renormalised
	    std::size_t index = s + res + r;
	    state[index].real *= factor;
	    state[index].imag *= factor;

	    // Get the indices that need to be zeroed out
	    index = s + not_res + r;
	    state[index].real = 0;
	    state[index].imag = 0;
	}
    }
    
}


template<typename Fp>
std::vector<typename Qubits<Type::Default, Fp>::Dist> Qubits<Type::Default, Fp>::generateDist()
{
    // Construct cumulative probability vector
    std::vector<Dist> dist;
    dist.push_back({0, 0});
    
    // Only stores non-zero probabilities
    for (std::size_t n = 0; n < dim; n++) {
	double prob = state[n].real * state[n].real
	    + state[n].imag * state[n].imag;
	if (prob != 0) {
	    /// \todo This line needs checking!
	    dist.push_back({n, dist.back().prob + prob});
	}
    }
    
    return dist;
}


template<typename Fp>
int Qubits<Type::Default, Fp>::measure(unsigned targ)
{
    // Calculate probabilty of measuring 0 on qubit targ
    double prob0 = prob(targ, 0);
    
    // Generate a random number and use it to generate the outcome
    // and calculate the renormalisation factor at the same time.
    // If the random number is less than the probability of getting
    // 0, then the outcome is 0, otherwise it is 1.
    double rand = random.getNum();
    int outcome;
    double factor;
    if (rand < prob0) {
	outcome = 0;
	factor = 1 / std::sqrt(prob0);
    }
    else {
	outcome = 1;
	factor = 1 / std::sqrt(1 - prob0);
    }
    
    // Collapse to outcome
    collapse(targ, outcome, factor);

    return outcome;
}


template<typename Fp>
std::size_t Qubits<Type::Default, Fp>::measureAll()
{
    // Construct cumulative probability vector
    std::vector<Dist> dist = generateDist();
    
    // Sample from the vector once
    std::size_t outcome = drawSample(dist);

    // Collapse to the outcome
    setBasisState(outcome);
    
    return outcome;
}


template<typename Fp>
double Qubits<Type::Default, Fp>::prob(unsigned targ, unsigned outcome)
    const
{
    double probability = 0;
    
    std::size_t k = 1 << targ;
    std::size_t res = outcome * k;
    
    for (std::size_t s = 0; s < dim; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the indices that need to be added
	    std::size_t index = s + res + r;

	    // Add the absolute value of the amplitude squared
	    probability += state[index].real * state[index].real
		+ state[index].imag * state[index].imag;
	}
    }

    return probability;
}

template<typename Fp>
double Qubits<Type::Default, Fp>::postselect(unsigned targ,
					     unsigned outcome)
{
    // Calculate renormalisation factor
    double probability = prob(targ, outcome);
    double factor = 1 / std::sqrt(probability);

    // Collapse to outcome
    collapse(targ, outcome, factor);

    return probability;
}


template<typename Fp>
std::size_t Qubits<Type::Default, Fp>::drawSample(const std::vector<Dist> & dist)
{
    std::size_t sampled = 0, L = 0, R = dist.size()-2, m = 0;
   
    // Generate a random number
    double rand = random.getNum();

    // Use a binary search to sample from the probability distribution
    while (L <= R) {
	m = (L+R)/2;
	if (rand <= dist[m].prob) {
	    R = m - 1;
	}
	else if (rand > dist[m+1].prob) {
	    L = m + 1;
	}
	else {
	    sampled = dist[m+1].index;
	    break;
	}
    }

    return sampled;
}


template<typename Fp>
std::vector<std::size_t> Qubits<Type::Default, Fp>::sample(unsigned targ,
							   std::size_t nsamples)
{
    // Initialise an empty vector, results[0] will store the amount of
    // times that 0 was sampled. Similarly for results[1]
    std::vector<std::size_t> results(2, 0);
    
    // Calculate probabilty of measuring 0 on qubit targ
    double prob0 = prob(targ, 0);
    
    // Generate a random number and use it to generate the outcome
    // If the random number is less than the probability of getting
    // 0, then the outcome is 0, otherwise it is 1.    
    for (std::size_t i = 0; i < nsamples; i++) {
	double rand = random.getNum();
	if (rand < prob0) {
	    results[0]++;
	}
	else {
	    results[1]++;
	}
    }

    return results;
}



template<typename Fp>
std::map<std::size_t, std::size_t> Qubits<Type::Default, Fp>::sampleAll(std::size_t nsamples)
{    
    // Construct cumulative probability vector
    std::vector<Dist> dist = generateDist();
    
    // Sample from the vector nmeas times
    std::map<std::size_t, std::size_t> results;
    for (std::size_t i = 0; i < nsamples; i++) {
	// The [] operator zero-initialises keys that don't exist.
	results[drawSample(dist)]++;
    }

    return results;
}


template<typename Fp>
std::map<std::size_t, std::size_t> Qubits<Type::Default, Fp>::sampleAll2(std::size_t nsamples)
{        
    // Generate a sorted list of nsamples random numbers
    std::vector<Fp> rand_list(nsamples);
    for (std::size_t i = 0; i < nsamples; i++) {
	rand_list[i] = random.getNum();
    }
    std::ranges::sort(rand_list);    
    //std::sort(rand_list.begin(), rand_list.end());
          
    std::map<std::size_t, std::size_t> results;
    
    Fp sum = 0;  // Running cumulative probability distribution
    std::size_t m = 0;  // Keeps track of rand_list index searched
    std::size_t count = 0;  // Keeps track of how many times an outcome has occured  
    
    /// \todo Can this be tidied up a bit and improved?
    for (std::size_t i = 0; i < dim; i++) {
	sum += state[i].real * state[i].real + state[i].imag * state[i].imag;
	count = 0;
	while (rand_list[m] < sum and m < nsamples) {
	    count++;
	    m++;
	}
	// Add to the map
	results[i] = count;
    }
        
    return results;
}

// Explicit instantiations
template class Qubits<Type::Default, float>;
template class Qubits<Type::Default, double>;
