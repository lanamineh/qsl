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
 * \file oneloop_test.cpp
 * \brief Speed test a single for loop (using masking)
 */

#include <vector>
#include <iostream>
#include <complex>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <random>
#include "qsl/utils/timer.hpp"
#include "qsl/utils/complex.hpp"
#include "qsl/utils/misc.hpp"

/**
 * \brief Apply the Pauli X gate to qubit number targ.
 */
void pauliX(std::vector<complex<double>> &state, std::uint8_t targ)
{
    std::size_t k = 1 << targ;
    // Create masks for breaking up the number
    std::size_t lower_mask = k - 1;
    std::size_t upper_mask = ~lower_mask;
    
    // Loop through all the other numbers
    for (std::size_t i = 0; i < state.size()/2; i++) {
	std::size_t lower = i & lower_mask;
	std::size_t upper = i & upper_mask;

	std::size_t index1 = lower + (upper << 1);
	std::size_t index2 = lower + k + (upper << 1);

	std::swap(state[index1], state[index2]);	
    }
}


/**
 * \brief Apply a phase shift to qubit number targ.
 */
void phaseShift(std::vector<complex<double>> &state, std::uint8_t targ, double angle)
{
    const complex<double> phase = complex<double>(std::cos(angle), std::sin(angle));
    
    std::size_t k = 1 << targ;
    // Create masks for breaking up the number
    std::size_t lower_mask = k - 1;
    std::size_t upper_mask = ~lower_mask;
    
    // Loop through all the other numbers
    for (std::size_t i = 0; i < state.size()/2; i++) {
	std::size_t lower = i & lower_mask;
	std::size_t upper = i & upper_mask;

	std::size_t index = lower + k + (upper << 1);

	complex temp = state[index];
	state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	state[index].imag = phase.real * temp.imag + phase.imag * temp.real;
    }
}

/**
 * \brief Rotate around the x-axis of the Bloch sphere. e^{-i*t*X/2}
 *
 * This single qubit gate applies the following 2x2 matrix to each
 * pair of |0> and |1> amplitudes for angle t:
 *
 * 
 *    -                        - 
 *   |   cos(t/2)   -i sin(t/2) |
 *   | -i sin(t/2)    cos(t/2)  |
 *    -                        - 
 *
 *
 */
void rotateX(std::vector<complex<double>> &state, std::uint8_t targ, double angle)
{
    // Store variables
    const double cos = std::cos(angle/2);
    const double sin = std::sin(angle/2);

    std::size_t k = 1 << targ;
    // Create masks for breaking up the number
    std::size_t lower_mask = k - 1;
    std::size_t upper_mask = ~lower_mask;
    
    // Loop through all the other numbers
    for (std::size_t i = 0; i < state.size()/2; i++) {
	std::size_t lower = i & lower_mask;
	std::size_t upper = i & upper_mask;

	std::size_t index_0 = lower + (upper << 1);
	std::size_t index_1 = lower + k + (upper << 1);
	
	// Store the values of |0> and |1> amplitudes
	complex<double> a0 = state[index_0];
	complex<double> a1 = state[index_1];
	
	// Write the new |0> amplitude
	state[index_0].real = a0.real * cos + a1.imag * sin;
	state[index_0].imag = a0.imag * cos - a1.real * sin;
	
	// Write the new |1> amplitude
	state[index_1].real = a1.real * cos + a0.imag * sin;
	state[index_1].imag = a1.imag * cos - a0.real * sin;
    }
}


/**
 * \brief Perform the CNOT gate on two qubits.
 */
void controlNot(std::vector<complex<double>> &state, std::uint8_t ctrl, std::uint8_t targ)
{
    // Find the bit positions of ctrl and targ
    std::size_t small_bit = 1 << std::min(ctrl, targ);
    std::size_t large_bit = 1 << std::max(ctrl, targ);

    // Create masks for the 3 sections that the bit string will be broken into 
    std::size_t lower_mask = small_bit - 1;
    std::size_t mid_mask = ((large_bit >> 1) - 1) ^ lower_mask;
    std::size_t upper_mask = ~(lower_mask | mid_mask);

    std::size_t targ_bit = (1 << targ);
    std::size_t ctrl_bit = (1 << ctrl);
    
    // Loop through all the other numbers with num_ones 1s using next
    // then break down that number into 3 parts to go on either side
    // of q1 and q2
    for (std::size_t i = 0; i < state.size()/4; i++) {
	std::size_t lower = i & lower_mask;
	std::size_t mid = i & mid_mask;
	std::size_t upper = i & upper_mask;

	// Calculate the index by adding together the 3 shifted sections
	// of x and small_bit and large_bit (which represent having 1
	// on ctrl and targ).
	std::size_t indexLo = lower + (mid << 1) + (upper << 2) + ctrl_bit;
	std::size_t indexUp = indexLo + targ_bit;

	std::swap(state[indexUp], state[indexLo]);
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
	    //pauliX(state_list[k], i);
	    //rotateX(state_list[k], i, phase_list[nqubits*k + i]);
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
