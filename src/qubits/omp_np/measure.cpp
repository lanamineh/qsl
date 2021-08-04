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
 * \file measure.cpp
 * \brief Contains the implementation of the measurement and sampling 
 *        for the number preserved Qubits class with OpenMP.
 */


#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <cmath>
#include <algorithm>


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::OmpNP, Fp>::collapse(unsigned targ, unsigned outcome,
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
#pragma omp parallel for num_threads(nthreads)
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
#pragma omp parallel for num_threads(nthreads)
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


template<std::floating_point Fp>
std::vector<typename qsl::Qubits<qsl::Type::OmpNP, Fp>::Dist> qsl::Qubits<qsl::Type::OmpNP, Fp>::generateDist()
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




template<std::floating_point Fp>
int qsl::Qubits<qsl::Type::OmpNP, Fp>::measure(unsigned targ)
{
    // Calculate probabilty of measuring 0 on qubit targ
    Fp prob0 = prob(targ, 0);
    
    // Generate a random number and use it to generate the outcome
    // and calculate the renormalisation factor at the same time.
    // If the random number is less than the probability of getting
    // 0, then the outcome is 0, otherwise it is 1.
    Fp rand = random.getNum();
    int outcome;
    Fp factor;
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


template<std::floating_point Fp>
std::size_t qsl::Qubits<qsl::Type::OmpNP, Fp>::measureAll()
{
    // Construct cumulative probability vector
    std::vector<Dist> dist = generateDist();
    
    // Sample from the vector once
    std::size_t outcome = drawSample(dist);

    // Collapse to the outcome
    setBasisState(outcome);
    
    return outcome;
}


template<std::floating_point Fp>
Fp qsl::Qubits<qsl::Type::OmpNP, Fp>::prob(unsigned targ, unsigned outcome) const
{
    Fp probability = 0;

    // Position and outcome of target qubit
    std::size_t k = 1 << targ;
    std::size_t res = outcome * k;

    // Create masks for breaking up the number
    std::size_t lower_mask = k - 1;
    std::size_t upper_mask = ~lower_mask;
        
    // Loop through all the other numbers
#pragma omp parallel for num_threads(nthreads) reduction(+:probability)
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


template<std::floating_point Fp>
Fp qsl::Qubits<qsl::Type::OmpNP, Fp>::postselect(unsigned targ,
						  unsigned outcome)
{
    // Calculate renormalisation factor
    Fp probability = prob(targ, outcome);
    Fp factor = 1 / std::sqrt(probability);

    // Collapse to outcome
    collapse(targ, outcome, factor);
    
    return probability;
}


template<std::floating_point Fp>
std::size_t qsl::Qubits<qsl::Type::OmpNP, Fp>::drawSample(const std::vector<Dist> & dist)
{
    std::size_t sampled = 0, L = 0, R = dist.size()-2, m = 0;
   
    // Generate a random number
    Fp rand = random.getNum();

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


template<std::floating_point Fp>
std::vector<std::size_t> qsl::Qubits<qsl::Type::OmpNP, Fp>::sample(unsigned targ,
								std::size_t nsamples)
{
    // Initialise an empty vector, results[0] will store the amount of
    // times that 0 was sampled. Similarly for results[1]
    //std::vector<std::size_t> results(2, 0);
    std::size_t num0 = 0; // number of measured zeros
    std::size_t num1 = 0; // number of measured ones
    
    // Calculate probabilty of measuring 0 on qubit targ
    Fp prob0 = prob(targ, 0);
    
    // Generate a random number and use it to generate the outcome
    // If the random number is less than the probability of getting
    // 0, then the outcome is 0, otherwise it is 1.    
#pragma omp parallel for num_threads(nthreads) reduction(+:num0,num1)
    for (std::size_t i = 0; i < nsamples; i++) {
	Fp rand = random.getNum();
	if (rand < prob0) {
	    //results[0]++;
	    num0++;
	}
	else {
	    //results[1]++;
	    num1++;
	}
    }

    return {num0,num1};
}


template<std::floating_point Fp>
std::map<std::size_t, std::size_t> qsl::Qubits<qsl::Type::OmpNP, Fp>::sampleAll(std::size_t nsamples)
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


template<std::floating_point Fp>
std::map<std::size_t, std::size_t> qsl::Qubits<qsl::Type::OmpNP, Fp>::sampleAll2(std::size_t nsamples)
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
template class qsl::Qubits<qsl::Type::OmpNP, float>;
template class qsl::Qubits<qsl::Type::OmpNP, double>;

