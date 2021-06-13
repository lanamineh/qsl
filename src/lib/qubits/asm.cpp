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
 * \file qubits_asm.cpp
 * \brief Contains the implementation of the Qubits class
 *
 */

#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <iostream>
#include <cmath>
#include <string>

const std::string Qubits<Type::Asm, double>::name = std::string("Qub<asm,double>");

Qubits<Type::Asm, double>::Qubits(unsigned nqubits_in) 
    : nqubits{ nqubits_in }, dim{ std::size_t(1) << nqubits },
      state{ dim }, random(0,1)

{
    // Make the all-zero state
    reset();   
}

Qubits<Type::Asm, double>::Qubits(const std::vector<complex<Fp_type>> & state)
    : nqubits{ checkStateSize(state) }, dim{ state.size() }, state{ state },
      random(0,1)
{
    //std::cout << "nqubits = " << nqubits << std::endl;
    //std::cout << "dim = " << dim << std::endl;
}

void Qubits<Type::Asm, double>::reset()
{
    for (std::size_t n=0; n<dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    state[0].real = 1;
}

void Qubits<Type::Asm, double>::setState(
    const std::vector<complex<Fp_type>> & state_in)
{
    if (state_in.size() != dim) {
	std::string msg = "Cannot assign state vector from different ";
	msg += "number of qubits";
	throw std::logic_error(msg);
    }

    // Set the new state vector
    state = state_in;
}

void Qubits<Type::Asm, double>::setBasisState(std::size_t index)
{
    // Clear the state - set all amplitudes to zero
    for (std::size_t n = 0; n < dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    // Set the amplitude for index to 1
    state[index].real = 1;
}



void Qubits<Type::Asm, double>::operator = (const Qubits & old) {

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
std::vector<complex<double>> Qubits<Type::Asm, double>::getState() const
{
    return state;
}

unsigned Qubits<Type::Asm, double>::getNumQubits() const
{
    return nqubits;
}


void Qubits<Type::Asm, double>::print() const
{
    std::cout << "Number of qubits = " << nqubits << std::endl;
    std::cout << state << std::endl;
}

extern "C" void
rotateX_x86(Qubits<Type::Asm> * ptr, unsigned targ, double angle);

extern "C" void
pauliX_x86(Qubits<Type::Asm> * ptr, unsigned targ);

void Qubits<Type::Asm, double>::rotateX(unsigned targ, double angle)
{
    rotateX_x86(this, targ, angle);
}

void Qubits<Type::Asm, double>::pauliX(unsigned targ)
{
    pauliX_x86(this, targ);
}

int Qubits<Type::Asm, double>::measure(unsigned targ)
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

    return outcome;
}

std::size_t Qubits<Type::Asm, double>::measureAll()
{
    // Construct cumulative probability vector
    std::vector<Dist> dist(dim + 1);
    dist[0].index = 0;
    dist[0].prob = 0;

    /// \todo Do not store zero probabilities
    for (std::size_t n = 0; n < dim; n++) {
	double prob = state[n].real * state[n].real
	    + state[n].imag * state[n].imag;
	dist[n + 1].index = n;
	dist[n + 1].prob = dist[n].prob + prob;
    }
    
    // Sample from the vector once
    std::size_t outcome = drawSample(dist);

    // Collapse to the outcome
    setBasisState(outcome);
    
    return outcome;
}

double Qubits<Type::Asm, double>::prob(unsigned targ, unsigned outcome)
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

double Qubits<Type::Asm, double>::postselect(unsigned targ, unsigned outcome)
{
    // Calculate renormalisation factor
    double probability = prob(targ, outcome);
    double factor = 1 / std::sqrt(probability);

    // Collapse to outcome and renormalise the state vectot simultaneously
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

    return probability;
}

std::size_t Qubits<Type::Asm, double>::drawSample(const std::vector<Dist> & dist)
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

std::vector<int> Qubits<Type::Asm, double>::sample(unsigned targ,
						   std::size_t nsamples)
{
    // Calculate probabilty of measuring 0 on qubit targ
    double prob0 = prob(targ, 0);
    
    // Generate a random number and use it to generate the outcome
    // and calculate the renormalisation factor at the same time.
    // If the random number is less than the probability of getting
    // 0, then the outcome is 0, otherwise it is 1.    
    std::vector<int> results(nsamples, 0);
    for (std::size_t i = 0; i < nsamples; i++) {
	double rand = random.getNum();
	if (rand > prob0) {
	    results[i] = 1;
	}
    }

    /// \todo Return the total numbers of 0s and 1s instead
    return results;
}

std::vector<std::size_t> Qubits<Type::Asm, double>::sampleAll(std::size_t nsamples)
{    
    // Construct cumulative probability vector
    std::vector<Dist> dist(dim + 1);
    dist[0].index = 0;
    dist[0].prob = 0;

    /// \todo Do not store zero probabilities
    for (std::size_t n = 0; n < dim; n++) {
	double prob = state[n].real * state[n].real
	    + state[n].imag * state[n].imag;
	dist[n + 1].index = n;
	dist[n + 1].prob = dist[n].prob + prob;
    }
    
    // Sample from the vector nmeas times
    std::vector<std::size_t> results(nsamples, 0);
    for (std::size_t i = 0; i < nsamples; i++) {
	results[i] = drawSample(dist);
    }

    /// \todo Return the counts of each sampled outcome instead
    return results;
}
