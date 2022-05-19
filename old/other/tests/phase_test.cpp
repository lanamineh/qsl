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

#include "cmake_defines.hpp"

/**
 * \brief Apply a phase shift to qubit number targ.
 */
void phaseShift(std::vector<qsl::complex<double>> &state, std::uint8_t targ, double angle)
{
    const qsl::complex<double> phase = qsl::complex<double>(std::cos(angle), std::sin(angle));
    
    std::size_t k = 1 << targ;
    for (std::size_t s = 0; s < state.size(); s += 2*k) { 
	for (std::size_t r = 0; r < k; r++) {
	    // Get the index of |1>
	    std::size_t index = s + k + r;

	    //state[index] *= phase;
	    qsl::complex<double> temp = state[index];
	    state[index].real = phase.real * temp.real - phase.imag * temp.imag;
	    state[index].imag = phase.real * temp.imag + phase.imag * temp.real;
	}
    }
}


int main()
{
    // Number of qubits and test length
    const std::uint8_t nqubits = NUM_QUBITS;
    const std::size_t test_length = TEST_LEN;
    
    std::cout << "Generating random vectors..." << std::endl;
    
    // Make a list of random state vectors
    std::vector<std::vector<qsl::complex<double>>> state_list;
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
