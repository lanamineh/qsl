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
 * \file test_flat.cpp
 * \brief Speed test for the flat array state vector
 *
 */

#include <vector>
#include <iostream>
#include <bitset>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <random>
#include "qsl/utils.hpp"

#include "cmake_defines.hpp"

/**
 * \brief Apply the Pauli X gate to qubit number targ.
 */
void pauliX(std::vector<double> &state, std::uint8_t targ)
{
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size()/2; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the indices that need to be switched
	    std::size_t index1 = 2*(s + r);
	    std::size_t index2 = 2*(s + k + r);

	    double tempr = state[index1];
	    double tempi = state[index1 + 1];
	    
	    state[index1] = state[index2];
	    state[index1 + 1] = state[index2 + 1];
	    
	    state[index2] = tempr;
	    state[index2 + 1] = tempi;
	}
    }
}


/**
 * \brief Apply a phase shift to qubit number targ.
 */
void phase(std::vector<double> &state, std::uint8_t targ, double angle)
{
    double phaser = std::cos(angle);
    double phasei = std::sin(angle);
    
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size()/2; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the index of |1>
	    std::size_t index = 2*(s + k + r);

	    //state[index] *= phase;
	    double tempr = state[index];
	    double tempi = state[index + 1];
	    
	    state[index] = phaser * tempr - phasei * tempi;
	    state[index + 1] = phaser * tempi + phasei * tempr;
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
void rotateX(std::vector<double> & state, std::uint8_t targ, double angle)
{
    // Store variables
    const double cos = std::cos(angle/2);
    const double sin = std::sin(angle/2);

    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size()/2; s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {

	    // Get the index of |0> and |1>
	    std::size_t index_0 = 2*(s + r);
	    std::size_t index_1 = 2*(s + k + r);

	    // Store the values of |0> and |1> amplitudes
	    double a0r = state[index_0];
	    double a0i = state[index_0+1];
	    double a1r = state[index_1];
	    double a1i = state[index_1+1];

	    // Write the new |0> amplitude
	    state[index_0] = a0r * cos + a1i * sin;
	    state[index_0+1] = a0i * cos - a1r * sin;

	    // Write the new |1> amplitude
	    state[index_1] = a1r * cos + a0i * sin;
	    state[index_1+1] = a1i * cos - a0r * sin;	    
	}
    }
}

/**
 * \brief Perform the CNOT gate on two qubits.
 */
void controlNot(std::vector<double> &state, std::uint8_t ctrl, std::uint8_t targ)
{
    std::size_t small_bit = 1 << std::min(ctrl, targ);
    std::size_t large_bit = 1 << std::max(ctrl, targ);

    std::size_t mid_incr = (small_bit << 1);
    std::size_t high_incr = (large_bit << 1);
    std::size_t targ_bit = (1 << targ);
    std::size_t ctrl_bit = (1 << ctrl);

    // Increment through the indices above largest bit (ctrl or targ)
    for(std::size_t i=0; i<state.size()/2; i+=high_incr) {
	// Increment through the middle set of bits
	for(std::size_t j=0; j<large_bit; j+=mid_incr) {
	    // Increment through the low set of bits
            for(std::size_t k=0; k<small_bit; k++) {
                // 2x2 matrix multiplication on the zero (i+j+k)
                // and one (i+j+k+targ_bit) indices. 
		//  mat_mul(op, state, i+j+k, i+j+k+targ_bit);
		std::size_t indexUp = 2*(i + j + k + ctrl_bit);
		std::size_t indexLo = indexUp + 2*targ_bit;

		double stateRealUp = state[indexUp];
                double stateImagUp = state[indexUp + 1];

                state[indexUp] = state[indexLo];
                state[indexUp + 1] = state[indexLo + 1];

                state[indexLo] = stateRealUp;
                state[indexLo + 1] = stateImagUp;

            }
	}
    }    
}

/**
 * \brief Normalise the state vector.
 */
double normalise(std::vector<double> &state)
{
    // Find the norm of the vector
    double norm = 0;
    for (std::size_t i = 0; i < state.size()/2; i++) {
	//norm += (state[i] * std::conj(state[i])).real();
	norm += state[2*i] * state[2*i] + state[2*i + 1] * state[2*i + 1];
    }
    norm = std::sqrt(norm);
   
    // Divide by the norm;
    for (std::size_t i = 0; i < state.size(); i++) {
	state[i] /= norm;
    }

    return norm;
}

/**
 * \brief Make a random state vector with nqubits
 */
std::vector<double> makeRandomState(std::uint8_t nqubits)
{
    std::size_t dim = 1 << nqubits;
    // Make the random int generator from -500 to 500
    std::random_device r;
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> distribution(-500,500);

    std::vector<double> state;
    for(std::size_t i=0; i<dim; i++) {
	double val_real = static_cast<double>(distribution(generator)) / 500;
	double val_imag = static_cast<double>(distribution(generator)) / 500;
	state.push_back(val_real);
	state.push_back(val_imag);
    }

    // Normalise the state vector
    normalise(state);
    
    return state;
}

int main()
{
    // Number of qubits and test length
    const std::uint8_t nqubits = NUM_QUBITS;
    const std::size_t test_length = TEST_LEN;

    std::cout << "Generating random vectors..." << std::endl;
    
    // Make a list of random state vectors
    // std::vector<std::vector<double>> state_list;
    // for(std::size_t k=0; k<test_length; k++) {
    // 	state_list.push_back(makeRandomState(nqubits));
    // }
    std::vector<double> state_list = makeRandomState(nqubits);

    // Make a list of random phases
    std::vector<double> phase_list{
	qsl::makeRandomPhases<double>(nqubits)
	    };

    std::cout << "Starting test..." << std::endl;
    qsl::Timer t;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	for(std::size_t i=0; i<nqubits; i++) {
#if GATE == 0
	    pauliX(state_list, i);
#elif GATE == 1
	    phase(state_list, i, phase_list[i]);
#elif GATE == 2
	    rotateX(state_list, i, phase_list[i]);
#endif
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
	    controlNot(state_list, i, i+1); 
	}
    }
    t.stop();
    std::cout << t.printElapsed() << std::endl;
    
    
    return 0;
}
