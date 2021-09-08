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
#include <complex>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <random>
#include "qsl/utils/timer.hpp"

/// SIMD type for one complex number (two packed 16-byte doubles)
using complex = double __attribute__((vector_size(32)));

/**
 * \brief Apply the Pauli X gate to qubit number targ.
 */
void pauliX(std::vector<complex> &state, std::uint8_t targ)
{
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the indices that need to be switched
	    std::size_t index1 = s + r;
	    std::size_t index2 = s + k + r;

	    complex temp = state[index1];
	    state[index1] = state[index2];
	    state[index2] = temp;
	}
    }
}

/**
 * \brief Apply a phase shift to qubit number targ.
 *
 * The state vector is stored as an array of 16-byte complex numbers, 
 * which can be stored in an xmm register as a pair of packed doubles, e.g.
 *
 *   [r, s].
 * 
 * The phase = a + ib can also be stored in two xmm registers as follows
 *
 *   [a, b] 
 * 
 * The result of multiplying the amplitude by the phase is
 *
 *   [v,w] = [a*r - b*s, a*s + b*r]
 *
 *  [ a*r, b*r]
 *  [-b*s, a*s]
 *  [  v ,  w ] + 
 *
 * which can be realised by 
 * 
 */
void phase(std::vector<complex> &state, std::uint8_t targ, double angle)
{
    //complex phase = complex(std::cos(angle), std::sin(angle));

    const complex phase_1{ std::cos(angle), std::cos(angle) }; 
    const complex phase_2{ -std::sin(angle), std::cos(angle) }; 
    
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the index of |1>
	    std::size_t index = s + k + r;

	    complex temp_r{ state[index][0], state[index][0] };
	    complex temp_s{ state[index][1], state[index][1] };
	    complex result_1{ phase_1 * temp_r };
	    complex result_2{ phase_2 * temp_s };
	    complex result{ result_1 + result_2 };
	    
	    state[index] = result;
	}
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
void rotateX(std::vector<complex> &state, std::uint8_t targ, double angle)
{
    // Store variables
    double cos = std::cos(angle/2);
    double sin = std::sin(angle/2);

    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {

	    // Get the index of |0> and |1>
	    std::size_t index_0 = s + r;
	    std::size_t index_1 = s + k + r;

	    // Store the values of |0> and |1> amplitudes
	    complex a0 = state[index_0];
	    complex a1 = state[index_1];

	    // Write the new |0> amplitude
	    state[index_0][0] = a0[0] * cos + a1[1] * sin;
	    state[index_0][1] = a0[1] * cos - a1[0] * sin;

	    // Write the new |1> amplitude
	    state[index_1][0] = a1[0] * cos + a0[1] * sin;
	    state[index_1][1] = a1[1] * cos - a0[0] * sin;
	    
	}
    }
}


/**
 * \brief Perform the CNOT gate on two qubits.
 */
void controlNot(std::vector<complex> &state, std::uint8_t ctrl, std::uint8_t targ)
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
		
		complex temp = state[indexUp];
		state[indexUp] = state[indexLo];
		state[indexLo] = temp;
            }
	}
    }    
}



/**
 * \brief Normalise the state vector.
 */
double normalise(std::vector<complex> &state)
{
    // Find the norm of the vector
    double norm = 0;
    for (std::size_t i = 0; i < state.size(); i++) {
	//norm += (state[i] * std::conj(state[i])).real();
	norm += (state[i][0] * state[i][0]) + (state[i][1] * state[i][1]);
    }
    norm = std::sqrt(norm);
   
    // Divide by the norm;
    for (std::size_t i = 0; i < state.size(); i++) {
	state[i][0] /= norm;
	state[i][1] /= norm;
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
std::vector<complex> makeRandomState(std::uint8_t nqubits)
{
    std::size_t dim = 1 << nqubits;
    // Make the random int generator from -500 to 500
    std::random_device r;
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> distribution(-500,500);

    std::vector<complex> state;
    for(std::size_t i=0; i<dim; i++) {
	double val_real = static_cast<double>(distribution(generator)) / 500;
	double val_imag = static_cast<double>(distribution(generator)) / 500;
	state.push_back(complex{val_real, val_imag});
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
    std::vector<std::vector<complex>> state_list;
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
	    rotateX(state_list[k], i, phase_list[nqubits*k + i]);
	    phase(state_list[k], i, phase_list[nqubits*k + i]);
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
