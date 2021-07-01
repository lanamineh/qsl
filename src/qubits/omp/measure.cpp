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
 * \file measure.cpp
 * \brief Contains the implementation of measurement and sampling
 *        in the OpenMP-based Qubits class.
 */

#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include <iostream>
#include <cmath>
#include <omp.h>

namespace qsl {

    template<std::floating_point Fp>
    void Qubits<Type::Omp, Fp>::collapse(unsigned targ, unsigned outcome,
					 Fp factor)
    {
	// Collapse to outcome and renormalise the state vector simultaneously   
	std::size_t k = 1 << targ;
	std::size_t res = outcome * k;
	std::size_t not_res = (outcome ^ 1) * k; 

	// Loop through the state vector and renormalise amplitudes associated with
	// the correct outcome and zero out those with the opposite outcome.
#pragma omp parallel for num_threads(nthreads)
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


    template<std::floating_point Fp>
    std::vector<typename Qubits<Type::Omp, Fp>::Dist> Qubits<Type::Omp, Fp>::generateDist()
    {
	// Construct cumulative probability vector
	std::vector<Dist> dist;
	dist.push_back({0, 0});
    
	// Only stores non-zero probabilities
	for (std::size_t n = 0; n < dim; n++) {
	    Fp prob = state[n].real * state[n].real
		+ state[n].imag * state[n].imag;
	    if (prob != 0) {
		/// \todo This line needs checking!
		dist.push_back({n, dist.back().prob + prob});
	    }
	}
    
	return dist;
    }


    template<std::floating_point Fp>
    int Qubits<Type::Omp, Fp>::measure(unsigned targ)
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


    template<std::floating_point Fp>
    std::size_t Qubits<Type::Omp, Fp>::measureAll()
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
    double Qubits<Type::Omp, Fp>::prob(unsigned targ, unsigned outcome)
	const
    {
	double probability = 0;
    
	std::size_t k = 1 << targ;
	std::size_t res = outcome * k;
    
#pragma omp parallel for num_threads(nthreads) 
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

    template<std::floating_point Fp>
    double Qubits<Type::Omp, Fp>::postselect(unsigned targ,
					     unsigned outcome)
    {
	// Calculate renormalisation factor
	double probability = prob(targ, outcome);
	double factor = 1 / std::sqrt(probability);

	// Collapse to outcome
	collapse(targ, outcome, factor);

	return probability;
    }


    template<std::floating_point Fp>
    std::size_t Qubits<Type::Omp, Fp>::drawSample(const std::vector<Dist> & dist)
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


    template<std::floating_point Fp>
    std::vector<std::size_t> Qubits<Type::Omp, Fp>::sample(unsigned targ,
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
#pragma omp parallel for num_threads(nthreads)
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



    template<std::floating_point Fp>
    std::map<std::size_t, std::size_t> Qubits<Type::Omp, Fp>::sampleAll(std::size_t nsamples)
    {    
	// Construct cumulative probability vector
	std::vector<Dist> dist = generateDist();
    
	// Sample from the vector nmeas times
	std::map<std::size_t, std::size_t> results;
#pragma omp parallel for num_threads(nthreads)
	for (std::size_t i = 0; i < nsamples; i++) {
	    // The [] operator zero-initialises keys that don't exist.
	    results[drawSample(dist)]++;
	}

	return results;
    }


    // Explicit instantiations
    template class Qubits<Type::Omp, float>;
    template class Qubits<Type::Omp, double>;

}
