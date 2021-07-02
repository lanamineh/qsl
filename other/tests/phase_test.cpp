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
 * \file phase_test.cpp
 * \brief Test the speed of the phase shift gate
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
#include "qsl/utils/complex.hpp"
#include "qsl/utils/misc.hpp"
#include "qsl/utils/quantum.hpp"
#include "qsl/utils/random.hpp"


/**
 * \brief Apply a phase shift to qubit number targ.
 */
void phaseShift(std::vector<complex<double>> &state, std::uint8_t targ, double angle)
{
    const complex<double> phase = complex<double>(std::cos(angle), std::sin(angle));
    
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the index of |1>
	    std::size_t index = s + k + r;

	    //state[index] *= phase;
	    complex<double> temp = state[index];
	    state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	    state[index].imag = phase.real * temp.imag + phase.imag * temp.real;
	}
    }
}


int main()
{
    std::uint8_t nqubits = 12;
    std::size_t dim = 1 << nqubits;
    
    std::cout << "Generating random vectors..." << std::endl;
    // Length of random tests
    std::size_t test_length = 20000;
    
    // Make a list of random state vectors
    std::vector<std::vector<complex<double>>> state_list;
    for(std::size_t k=0; k<test_length; k++) {
	state_list.push_back(makeRandomState(nqubits));
    }

    // Make a list of random phases
    std::vector<double> phase_list = makeRandomPhases(test_length*nqubits);

    // QuEST test ----------------------------------------------------
    QuESTEnv env = createQuESTEnv();

    // Make a list of random state vectors
    std::vector<Qureg> qureg_list;
    for(std::size_t k=0; k<test_length; k++) {
	Qureg qubits = createQureg(nqubits, env);
	for (std::size_t i = 0; i < dim; i++) {
	    qubits.stateVec.real[i] = state_list[k][i].real;
	    qubits.stateVec.imag[i] = state_list[k][i].imag;
	}
	qureg_list.push_back(qubits);
    }
    
    std::cout << "Starting test..." << std::endl;
    qsl::Timer t;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	//for(int i=0; i<nqubits; i++) {
	//    phaseShift(state_list[k], i, phase_list[nqubits*k + i]);
	//}
	phaseShift(state_list[k], 0, phase_list[k]);
    }
    t.stop();
    std::cout << t.printElapsed() << std::endl;

    std::cout << "Starting quest test..." << std::endl;
    t.start();
    // Run a speed test for all the state vectors
    for(std::size_t k=0; k<test_length; k++) {
	// Apply Pauli X and phase shift to all qubits
	//for(std::size_t i=0; i<nqubits; i++) {
	//   phaseShift(qureg_list[k], i, phase_list[nqubits*k + i]);
	//}
	phaseShift(qureg_list[k], 0, phase_list[k]);
    }
    t.stop();
    std::cout << t.printElapsed() << std::endl;


    for(std::size_t k=0; k<test_length; k++) {
	destroyQureg(qureg_list[k], env);
    }
    destroyQuESTEnv(env);
    
    
    return 0;
}