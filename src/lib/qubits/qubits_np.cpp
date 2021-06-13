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
 * \file qubits_np.cpp
 * \brief Contains the implementation of the number preserved Qubits class.
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
const std::string Qubits<Type::NP, double>::name =
    std::string("Qub<np,double>");

template<>
const std::string Qubits<Type::NP, float>::name =
    std::string("Qub<np,float>");

template<typename Fp>
Qubits<Type::NP, Fp>::Qubits(unsigned nqubits_in, unsigned nones_in) 
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


/**
 * \brief Check if the state is number preserving and find the number
 * of ones if it is.
 *
 * This function is probably inefficient - needs improving.
 * Also might make it a member function of the class 
 */
template<typename Fp>
unsigned checkStateNP(const std::vector<complex<Fp>> & state)
{
    unsigned nones = 0;
    bool found = false;
    
    for (std::size_t i = 0; i < state.size(); i++) {
	Fp amp = state[i].real * state[i].real +
	    state[i].imag * state[i].imag;
	// If amplitude is non-zero store the number of ones
	if (amp != 0) {
	    unsigned weight = hammingWeight(i);
	    if (found == false) {
		nones = weight;
		found = true;
	    }
	    else if (nones != weight) {
		throw std::logic_error("Input state is not number preserving.");
	    }
	}
    }

    return nones;
}



template<typename Fp>
Qubits<Type::NP, Fp>::Qubits(const std::vector<complex<Fp>> & state)
    : nqubits{ checkStateSize(state) }, dim{ state.size() },
      nones{ checkStateNP(state) }, state{ state }, random(0,1)
{
    if (nones < 1 or nones >= nqubits) {
	throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
    }

    initLookup();
    //std::cout << "nqubits = " << nqubits << std::endl;
    //std::cout << "dim = " << dim << std::endl;
}


template<typename Fp>
void Qubits<Type::NP, Fp>::reset()
{
    for (std::size_t n = 0; n < dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    
    // Starting number (nones 1s in the least significant positions)
    std::size_t idx = (1ULL << nones) - 1ULL;
    state[idx].real = 1;
}

template<typename Fp>
void Qubits<Type::NP, Fp>::setState(const std::vector<complex<Fp>> & state_in)
{   
    if (state_in.size() != dim) {
	std::string msg = "Cannot assign state vector from different ";
	msg += "number of qubits";
	throw std::logic_error(msg);
    }

    nones = checkStateNP(state_in);
    if (nones < 1 or nones >= nqubits) {
	throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
    }

    initLookup();
    
    // Set the new state vector
    state = state_in;
}

template<typename Fp>
void Qubits<Type::NP, Fp>::setBasisState(std::size_t index)
{
    unsigned nones_old = nones;
    nones = hammingWeight(index);
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



template<typename Fp>
void Qubits<Type::NP, Fp>::operator = (const Qubits & old)
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
std::vector<complex<Fp>> Qubits<Type::NP, Fp>::getState() const
{
    return state;
}

template<typename Fp>
unsigned Qubits<Type::NP, Fp>::getNumQubits() const
{
    return nqubits;
}


template<typename Fp>
unsigned Qubits<Type::NP, Fp>::getNumOnes() const
{
    return nones;
}


template<typename Fp>
void Qubits<Type::NP, Fp>::setNumOnes(unsigned nones_in) 
{
    if (nones_in < 1 or nones_in >= nqubits) {
	throw std::logic_error("The number of ones must be 1 <= nones < nqubits.");
    }
    
    nones = nones_in;
    initLookup();
    reset();
}


/// Generate the lookup vector for a given bit-string length and number of ones
[[nodiscard]] std::vector<std::size_t> generateLookup(unsigned len, unsigned ones)
{
    std::vector<std::size_t> list;
    // Starting number (ones in the least significant positions)
    std::size_t x = (1ULL << ones) - 1ULL;
    // Last number (ones in the most significant positions)
    std::size_t end = x << (len - ones);

    // Loop through all the other numbers
    while(x <= end) {
	list.push_back(x);
	next(x);
    }

    return list;
}


template<typename Fp>
void Qubits<Type::NP, Fp>::initLookup() 
{
    lookup.clear();
    // Generate lookup tables
    // For use in 1-qubit gates and wavefunction collapse
    lookup[{nqubits-1, nones-1}] = generateLookup(nqubits-1, nones-1);
    // Used for 2-qubit gates
    lookup[{nqubits-2, nones-1}] = generateLookup(nqubits-2, nones-1);
    lookup[{nqubits-2, nones-2}] = generateLookup(nqubits-2, nones-2);
    // Used in measurement functions - e.g. generate distribution, collapse
    lookup[{nqubits-1, nones}] = generateLookup(nqubits-1, nones);
    lookup[{nqubits, nones}] = generateLookup(nqubits, nones);   
}



template<typename Fp>
void Qubits<Type::NP, Fp>::print() const
{
    std::cout << "Number of qubits = " << nqubits << std::endl;
    std::cout << "Number of ones = " << nones << std::endl;
    std::cout << state << std::endl;
}


template<typename Fp>
void Qubits<Type::NP, Fp>::phase(unsigned targ, Fp angle)
{
    complex<Fp> phase = complex(std::cos(angle), std::sin(angle));

    // Position of target qubit
    std::size_t k = 1 << targ;

    // Create masks for breaking up the number
    std::size_t lower_mask = k - 1;
    std::size_t upper_mask = ~lower_mask;
    
    // Loop through all the other numbers
//#pragma omp parallel for
    for (std::size_t i = 0; i < lookup.at({nqubits-1, nones-1}).size(); i++) {
	std::size_t x = lookup.at({nqubits-1, nones-1})[i];
	std::size_t lower = x & lower_mask;
	std::size_t upper = x & upper_mask;
	// Get index for 1 in target position
	std::size_t index = lower + k + (upper << 1);

	complex temp = state[index];
	state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	state[index].imag = phase.real * temp.imag + phase.imag * temp.real;
    }
}


template<typename Fp>
void Qubits<Type::NP, Fp>::controlPhase(unsigned ctrl, unsigned targ,
					Fp angle)
{
    // Gate does nothing if there is only one 1.
    if (nones < 2) {
	return;
    }
    
    complex<Fp> phase = complex(std::cos(angle), std::sin(angle));

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
    for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-2}).size(); i++) {
	std::size_t x = lookup.at({nqubits-2, nones-2})[i];
	std::size_t lower = x & lower_mask;
	std::size_t mid = x & mid_mask;
	std::size_t upper = x & upper_mask;

	// Calculate the index by adding together the 3 shifted sections
	// of x and small_bit and large_bit (which represent having 1
	// on ctrl and targ).
	std::size_t index = lower + small_bit + (mid << 1) + large_bit + (upper << 2);

	complex amp = state[index];
	state[index].real = phase.real * amp.real - phase.imag * amp.imag;
	state[index].imag = phase.real * amp.imag + phase.imag * amp.real;
    }
}


template<typename Fp>
void Qubits<Type::NP, Fp>::swap(unsigned q1, unsigned q2)
{
    // Find the bit positions of ctrl and targ
    std::size_t small_bit = 1 << std::min(q1, q2);
    std::size_t large_bit = 1 << std::max(q1, q2);

    // Create masks for the 3 sections that the bit string will be broken into 
    std::size_t lower_mask = small_bit - 1;
    std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
    std::size_t upper_mask = ~(lower_mask | mid_mask);
        
    // Loop through all the other numbers with num_ones 1s 
    // then break down that number into 3 parts to go on either side
    // of q1 and q2
    for (std::size_t i = 0; i < lookup.at({nqubits-2, nones-1}).size(); i++) {
	std::size_t x = lookup.at({nqubits-2, nones-1})[i];
	std::size_t lower = x & lower_mask;
	std::size_t mid = x & mid_mask;
	std::size_t upper = x & upper_mask;

	// Calculate the index by adding together the 3 shifted sections
	// of x and small_bit and large_bit (which represent having 1
	// on ctrl and targ).
	std::size_t index01 = lower + small_bit + (mid << 1) + (upper << 2);
	std::size_t index10 = lower + (mid << 1) + large_bit + (upper << 2);

	std::swap(state[index01], state[index10]);
    }
}


template<typename Fp>
void Qubits<Type::NP, Fp>::collapse(unsigned targ, unsigned outcome,
				    Fp factor)
{
    // Collapse to outcome and renormalise the state vector simultaneously   
    std::size_t k = 1 << targ;
    std::size_t res = outcome * k;
    std::size_t not_res = (outcome ^ 1) * k; 

    // Create masks for breaking up the number
    std::size_t lower_mask = k - 1;
    std::size_t upper_mask = ~lower_mask;
    
    // Loop through state vector to renormalise amplitides associated
    // with the correct outcome
    for (std::size_t i = 0; i < lookup.at({nqubits-1, nones-outcome}).size(); i++) {
	std::size_t x = lookup.at({nqubits-1, nones-outcome})[i];
	std::size_t lower = x & lower_mask;
	std::size_t upper = x & upper_mask;
	// Get index for outcome in target position
	std::size_t index = lower + res + (upper << 1);

	state[index].real *= factor;
	state[index].imag *= factor;
    }

    // Now zero out the remaining amplitudes corresponding to the opposite outcome
    for (std::size_t i = 0; i < lookup.at({nqubits-1, nones-(outcome^1)}).size(); i++) {
	std::size_t x = lookup.at({nqubits-1, nones-(outcome^1)})[i];
	std::size_t lower = x & lower_mask;
	std::size_t upper = x & upper_mask;
	// Get index for outcome in target position
	std::size_t index = lower + not_res + (upper << 1);

	state[index].real = 0;
	state[index].imag = 0;
    }
    
}


template<typename Fp>
std::vector<typename Qubits<Type::NP, Fp>::Dist> Qubits<Type::NP, Fp>::generateDist()
{
    // Construct cumulative probability vector
    std::vector<Dist> dist;
    dist.push_back({0, 0});

    for (std::size_t i = 0; i < lookup.at({nqubits, nones}).size(); i++) {
	std::size_t x = lookup.at({nqubits, nones})[i];
	Fp prob = state[x].real * state[x].real
	    + state[x].imag * state[x].imag;
	/// \todo This line needs checking!
	if (prob != 0) {
	    dist.push_back({x, dist.back().prob + prob});
	}
    }

    return dist;
}




template<typename Fp>
int Qubits<Type::NP, Fp>::measure(unsigned targ)
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
std::size_t Qubits<Type::NP, Fp>::measureAll()
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
double Qubits<Type::NP, Fp>::prob(unsigned targ, unsigned outcome) const
{
    double probability = 0;

    // Position and outcome of target qubit
    std::size_t k = 1 << targ;
    std::size_t res = outcome * k;

    // Create masks for breaking up the number
    std::size_t lower_mask = k - 1;
    std::size_t upper_mask = ~lower_mask;
        
    // Loop through all the other numbers
    for (std::size_t i = 0; i < lookup.at({nqubits-1, nones-outcome}).size(); i++) {
	std::size_t x = lookup.at({nqubits-1, nones-outcome})[i];
	std::size_t lower = x & lower_mask;
	std::size_t upper = x & upper_mask;
	// Get index for outcome in target position
	std::size_t index = lower + res + (upper << 1);

	// Add the absolute value of the amplitude squared
	probability += state[index].real * state[index].real
	    + state[index].imag * state[index].imag;
    }

    return probability;
}


template<typename Fp>
double Qubits<Type::NP, Fp>::postselect(unsigned targ,
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
std::size_t Qubits<Type::NP, Fp>::drawSample(const std::vector<Dist> & dist)
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
std::vector<std::size_t> Qubits<Type::NP, Fp>::sample(unsigned targ,
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
std::map<std::size_t, std::size_t> Qubits<Type::NP, Fp>::sampleAll(std::size_t nsamples)
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
std::map<std::size_t, std::size_t> Qubits<Type::NP, Fp>::sampleAll2(std::size_t nsamples)
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

    for (std::size_t i = 0; i < lookup.at({nqubits, nones}).size(); i++) {
	std::size_t x = lookup.at({nqubits, nones})[i];
	sum += state[x].real * state[x].real + state[x].imag * state[x].imag;
	count = 0;
	while (rand_list[m] < sum and m < nsamples) {
	    count++;
	    m++;
	}
	// Add to the map
	results[x] = count;
    }

    return results;
}

// Explicit instantiations
template class Qubits<Type::NP, float>;
template class Qubits<Type::NP, double>;
