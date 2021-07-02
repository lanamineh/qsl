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
#include "qsl/utils/timer.hpp"
#include "qsl/utils/misc.hpp"

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
    std::uint8_t nqubits = 12;
    //std::size_t dim = 1 << nqubits;
    
    std::cout << "Generating random vectors..." << std::endl;
    // Length of random tests
    std::size_t test_length = 20000;
    
    // Make a list of random state vectors
    std::vector<State> state_list;
    for(std::size_t k=0; k<test_length; k++) {
	state_list.push_back(makeRandomState(nqubits));
    }

    // Make a list of random phases
    std::vector<double> phase_list;
    for(std::size_t k=0; k<test_length*nqubits; k++) {
	phase_list.push_back(makeRandomNumber(-M_PI, M_PI));
    }

    std::cout << "Starting 1-qubit gate test..." << std::endl;
    qsl::Timer t;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	for(std::size_t i=0; i<nqubits; i++) {
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