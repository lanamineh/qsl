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
 * \file test_separate.cpp
 *
 * \brief We construct the state vector using two separate vectors for real
 * and imaginary components in a struct. This is how quest do their state
 * vectors.
 */

#include <vector>
#include <iostream>
#include <complex>
#include <bitset>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <random>
#include "qsl/utils.hpp"

#include "cmake_defines.hpp"

struct State
{
    std::vector<double> real;
    std::vector<double> imag;
};

/**
 * \brief Apply the Pauli X gate to qubit number targ.
 */
void pauliX(State &state, std::uint8_t targ)
{  
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.real.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the indices that need to be switched
	    std::size_t index1 = (s + r);
	    std::size_t index2 = (s + k + r);

	    double tempr = state.real[index1];
	    double tempi = state.imag[index1];
	    
	    state.real[index1] = state.real[index2];
	    state.imag[index1] = state.imag[index2];
	    
	    state.real[index2] = tempr;
	    state.imag[index2] = tempi;
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
void rotateX(State & state, std::uint8_t targ, double angle)
{
    // Store variables
    const double cos = std::cos(angle/2);
    const double sin = std::sin(angle/2);

    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.real.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {

	    // Get the index of |0> and |1>
	    std::size_t index_0 = s + r;
	    std::size_t index_1 = s + k + r;

	    // Store the values of |0> and |1> amplitudes
	    double a0r = state.real[index_0];
	    double a0i = state.imag[index_0];
	    double a1r = state.real[index_1];
	    double a1i = state.imag[index_1];

	    // Write the new |0> amplitude
	    state.real[index_0] = a0r * cos + a1i * sin;
	    state.imag[index_0] = a0i * cos - a1r * sin;

	    // Write the new |1> amplitude
	    state.real[index_1] = a1r * cos + a0i * sin;
	    state.imag[index_1] = a1i * cos - a0r * sin;	    
	}
    }
}


/**
 * \brief Apply a phase shift to qubit number targ.
 */
void phaseShift(State &state, std::uint8_t targ, double angle)
{
    double phaser = std::cos(angle);
    double phasei = std::sin(angle);
    
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.real.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the index of |1>
	    std::size_t index = s + k + r;

	    //state[index] *= phase;
	    double tempr = state.real[index];
	    double tempi = state.imag[index];
	    
	    state.real[index] = phaser * tempr - phasei * tempi;
	    state.imag[index] = phaser * tempi + phasei * tempr;
	}
    }
}


/**
 * \brief Perform the CNOT gate on two qubits.
 */
void controlNot(State &state, std::uint8_t ctrl, std::uint8_t targ)
{
    std::size_t small_bit = 1 << std::min(ctrl, targ);
    std::size_t large_bit = 1 << std::max(ctrl, targ);

    std::size_t mid_incr = (small_bit << 1);
    std::size_t high_incr = (large_bit << 1);
    std::size_t targ_bit = (1 << targ);
    std::size_t ctrl_bit = (1 << ctrl);

    // Increment through the indices above largest bit (ctrl or targ)
    for(std::size_t i=0; i<state.real.size(); i+=high_incr) {
	// Increment through the middle set of bits
	for(std::size_t j=0; j<large_bit; j+=mid_incr) {
	    // Increment through the low set of bits
            for(std::size_t k=0; k<small_bit; k++) {
                // 2x2 matrix multiplication on the zero (i+j+k)
                // and one (i+j+k+targ_bit) indices. 
		//  mat_mul(op, state, i+j+k, i+j+k+targ_bit);
		std::size_t indexUp = i + j + k + ctrl_bit;
		std::size_t indexLo = indexUp + targ_bit;
		
		double stateRealUp = state.real[indexUp];
                double stateImagUp = state.real[indexUp];

                state.real[indexUp] = state.real[indexLo];
                state.imag[indexUp] = state.imag[indexLo];

                state.real[indexLo] = stateRealUp;
                state.imag[indexLo] = stateImagUp;
            }
	}
    }
    
}


/**
 * \brief Normalise the state vector.
 */
double normalise(State &state)
{
    // Find the norm of the vector
    double norm = 0;
    for (std::size_t i = 0; i < state.real.size(); i++) {
	//norm += (state[i] * std::conj(state[i])).real();
	norm += state.real[i] * state.real[i] + state.imag[i] * state.imag[i];
    }
    norm = std::sqrt(norm);
   
    // Divide by the norm;
    for (std::size_t i = 0; i < state.real.size(); i++) {
	state.real[i] /= norm;
	state.imag[i] /= norm;
    }

    return norm;
}

/**
 * \brief Make a random state vector with nqubits
 */
State makeRandomState(std::uint8_t nqubits)
{
    std::size_t dim = 1 << nqubits;
    // Make the random int generator from -500 to 500
    std::random_device r;
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> distribution(-500,500);

    State state;
    for(std::size_t i=0; i<dim; i++) {
	double val_real = static_cast<double>(distribution(generator)) / 500;
	double val_imag = static_cast<double>(distribution(generator)) / 500;
	state.real.push_back(val_real);
	state.imag.push_back(val_imag);
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
    std::vector<State> state_list;
    for(std::size_t k=0; k<test_length; k++) {
	state_list.push_back(makeRandomState(nqubits));
    }

    // Make a list of random phases
    std::vector<double> phase_list{
	qsl::makeRandomPhases<double>(test_length * nqubits)
	    };
    
    std::cout << "Starting 1-qubit gate test..." << std::endl;
    qsl::Timer t;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	for(std::size_t i=0; i<nqubits; i++) {
	    //pauliX(state_list[k], i);
	    //phaseShift(state_list[k], i, phase_list[nqubits*k + i]);
	    rotateX(state_list[k], i, phase_list[nqubits*k + i]);
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
