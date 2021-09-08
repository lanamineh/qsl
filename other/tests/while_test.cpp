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
 * \file test_struct.cpp
 * \brief Speed test our own complex number struct
 *
 */

#include <vector>
#include <iostream>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <random>
#include "qsl/utils/timer.hpp"
#include "qsl/utils/complex.hpp"
#include "qsl/utils/misc.hpp"

enum {
    READ = 0,
    WRITE = 1,
};

/**
 * \brief Apply the Pauli X gate to qubit number targ.
 */
void pauliX(std::vector<complex<double>> &state, std::uint8_t targ)
{
    std::size_t k = 1 << targ;

    std::size_t s_step = 2*k;
    std::size_t s_limit = state.size();
    std::size_t s = 0;
    while(s < s_limit) {
	
	std::size_t r = 0;
	while(r < k) {
	    
	    // Get the indices that need to be switched
	    std::size_t index1 = s + r;
	    std::size_t index2 = s + k + r++; // Note r post-increment

	    complex<double> temp = state[index1];
	    state[index1] = state[index2];
	    state[index2] = temp;

	}

	s += s_step;
    }
}


/**
 * \brief Apply a phase shift to qubit number targ.
 */
void phaseShift(std::vector<complex<double>> &state, std::uint8_t targ, double angle)
{
    complex<double> phase = complex<double>(std::cos(angle), std::sin(angle));
    
    std::size_t k = 1 << targ;

    std::size_t s_step = 2*k;
    std::size_t s_limit = state.size();
    std::size_t s = 0;
    while(s < s_limit) {

	std::size_t prefetch = s+s_step;
	__builtin_prefetch(&state[prefetch]);	
	
	std::size_t r = 0;
	while(r < k) {
	    
	    // Get the index of |1>
	    std::size_t index = s + k + r++; // Note r post-increment

	    //state[index] *= phase;
	    complex<double> temp = state[index];
	    state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	    state[index].imag = phase.real * temp.imag + phase.imag * temp.real;

	}

	__builtin_prefetch(&state[prefetch]);	
	
	s += s_step;
    }
}


/**
 * \brief Perform the CNOT gate on two qubits.
 */
void controlNot(std::vector<complex<double>> &state, std::uint8_t ctrl, std::uint8_t targ)
{
    std::size_t small_bit = 1 << std::min(ctrl, targ);
    std::size_t large_bit = 1 << std::max(ctrl, targ);

    std::size_t mid_incr = (small_bit << 1);
    std::size_t high_incr = (large_bit << 1);
    std::size_t targ_bit = (1 << targ);
    std::size_t ctrl_bit = (1 << ctrl);

    // Increment through the indices above largest bit (ctrl or targ)
    for(std::size_t i=0; i<state.size(); i+=high_incr) {
	// Increment through the middle set of bits
	for(std::size_t j=0; j<large_bit; j+=mid_incr) {
	    // Increment through the low set of bits
            for(std::size_t k=0; k<small_bit; k++) {
                // 2x2 matrix multiplication on the zero (i+j+k)
                // and one (i+j+k+targ_bit) indices. 
		//  mat_mul(op, state, i+j+k, i+j+k+targ_bit);
		std::size_t indexUp = i + j + k + ctrl_bit;
		std::size_t indexLo = indexUp + targ_bit;
		
		complex<double> temp = state[indexUp];
		state[indexUp] = state[indexLo];
		state[indexLo] = temp;
            }
	}
    }    
}



/**
 * \brief Normalise the state vector.
 */
double normalise(std::vector<complex<double>> &state)
{
    // Find the norm of the vector
    double norm = 0;
    for (std::size_t i = 0; i < state.size(); i++) {
	//norm += (state[i] * std::conj(state[i])).real();
	norm += state[i].real * state[i].real + state[i].imag * state[i].imag;
    }
    norm = std::sqrt(norm);
   
    // Divide by the norm;
    for (std::size_t i = 0; i < state.size(); i++) {
	state[i].real /= norm;
	state[i].imag /= norm;
    }

    return norm;
}

/**
 * \brief Generate a random number between a and b
 */
double makeRandomNumber(double a, double b) {

    // Make the random int generator from -500 to 500
    std::random_device r;
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> distribution(-500,500);

    double val = static_cast<double>(distribution(generator));
    
    // Return the number
    double result =  (b+a)/2 + val*(b-a)/1000;

    return result;
}

/**
 * \brief Make a random state vector with nqubits
 */
std::vector<complex<double>> makeRandomState(std::uint8_t nqubits)
{
    std::size_t dim = 1 << nqubits;
    // Make the random int generator from -500 to 500
    std::random_device r;
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> distribution(-500,500);

    std::vector<complex<double>> state;
    for(std::size_t i=0; i<dim; i++) {
	double val_real = static_cast<double>(distribution(generator)) / 500;
	double val_imag = static_cast<double>(distribution(generator)) / 500;
	state.push_back(complex<double>(val_real, val_imag));
    }

    // Normalise the state vector
    normalise(state);
    
    return state;
}


int main()
{
    std::uint8_t nqubits = 12;
    //std::size_t dim = 1 << nqubits;
    
    std::cout << "Generating random vectors..." << std::endl;
    // Length of random tests
    std::size_t test_length = 20000;
    
    // Make a list of random state vectors
    std::vector<std::vector<complex<double>>> state_list;
    for(std::size_t k=0; k<test_length; k++) {
	state_list.push_back(makeRandomState(nqubits));
    }

    // Make a list of random phases
    std::vector<double> phase_list;
    for(std::size_t k=0; k<test_length*nqubits; k++) {
	phase_list.push_back(makeRandomNumber(-M_PI, M_PI));
    }

    std::cout << "Starting test..." << std::endl;
    qsl::Timer t;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	for(int i=0; i<nqubits; i++) {
	    pauliX(state_list[k], i);
	    phaseShift(state_list[k], i, phase_list[nqubits*k + i]);
	}
    }
    t.stop();
    std::cout << t.printElapsed() << std::endl;

    std::cout << "Starting 2-qubit gate test..." << std::endl;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	for(int i=0; i<nqubits-1; i++) {
	    controlNot(state_list[k], i, i+1); 
	}
    }
    t.stop();
    std::cout << t.printElapsed() << std::endl;
    

    
    return 0;
}
