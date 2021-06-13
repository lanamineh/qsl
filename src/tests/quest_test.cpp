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
 * \file quest_test.cpp
 * \brief Testing std::complex with Quest
 *
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
#include <QuEST.h>

using complex = std::complex<double>;

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
 */
void phaseShift(std::vector<complex> &state, std::uint8_t targ, double angle)
{
    complex phase = complex(std::cos(angle), std::sin(angle));
    
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the index of |1>
	    std::size_t index = s + k + r;

	    state[index] *= phase;
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
	norm += (state[i] * std::conj(state[i])).real();
    }
    norm = std::sqrt(norm);
   
    // Divide by the norm;
    for (std::size_t i = 0; i < state.size(); i++) {
	state[i] /= norm;
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
	state.push_back(complex(val_real, val_imag));
    }

    // Normalise the state vector
    normalise(state);
    
    return state;
}


int main()
{
    std::uint8_t nqubits = 12;
    std::size_t dim = 1 << nqubits;

    
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
	for(std::size_t i=0; i<nqubits; i++) {
	    pauliX(state_list[k], i);
	    phaseShift(state_list[k], i, phase_list[nqubits*k + i]);
	}
    }
    t.stop();
    std::cout << t.printElapsed() << std::endl;



    // QuEST test ----------------------------------------------------
    QuESTEnv env = createQuESTEnv();

    // Make a list of random state vectors
    std::vector<Qureg> qureg_list;
    for(std::size_t k=0; k<test_length; k++) {
	Qureg qubits = createQureg(nqubits, env);
	for (std::size_t i = 0; i < dim; i++) {
	    qubits.stateVec.real[i] = state_list[k][i].real();
	    qubits.stateVec.imag[i] = state_list[k][i].imag();
	}
	qureg_list.push_back(qubits);
    }


    std::cout << "Starting quest test..." << std::endl;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	for(std::size_t i=0; i<nqubits; i++) {
	    //pauliX(qureg_list[k], i);
	    phaseShift(qureg_list[k], i, phase_list[nqubits*k + i]);
	}
    }
    t.stop();
    std::cout << t.printElapsed() << std::endl;


    for(std::size_t k=0; k<test_length; k++) {
	destroyQureg(qureg_list[k], env);
    }
    destroyQuESTEnv(env);
    
    
    return 0;
}
