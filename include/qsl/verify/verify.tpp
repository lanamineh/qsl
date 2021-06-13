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
 * \file verify.cpp
 * \brief Implementation of gate verification against quest
 *
 *
 */

#include "qsl/qubits.hpp"
#include "qsl/quest.hpp"
#include "qsl/utils/quantum.hpp"
#include "qsl/utils/misc.hpp"
#include "qsl/utils/timer.hpp"
#include <iostream>
#include <functional>
#include <map>

///\todo Fix this for float
template<Type T, typename Fp = double>
void verifyPostselect(unsigned nqubits)
{
    std::cout << "============================================" << std::endl;
    std::cout << "Verifying postselect on " << nqubits << " qubits."
	      << std::endl;

    std::vector<complex<Fp>> state = makeRandomState(nqubits);
    Qubits<T, Fp> qubits(state);
    Quest quest(state);

    std::cout << "Printing probabilities of collapsing to outcome" << std::endl;
    
    for (unsigned n = 0; n < nqubits; n++) {

	// Make local copy
	Qubits<T> q1(qubits); // Make a local copy
	Quest q2(quest); // Make a local copy
	
	// Postselect on 0
	double qubits_0 = q1.postselect(n, 0);
	double quest_0 = q2.postselect(n, 0);
	double norm_0 = norm(q1.getState());
	double distance_0 = fubiniStudy(q1.getState(), q2.getState());

	// Reset the state vectors
	q1 = qubits;
	q2 = quest;

	// Postselect on 1
	double qubits_1 = q1.postselect(n, 1);
	double quest_1 = q2.postselect(n, 1);
	double norm_1 = norm(q1.getState());
	double distance_1 = fubiniStudy(q1.getState(), q2.getState());
	
	std::cout << "Qubit " << n
		  << ": Qubits = [" << qubits_0 << ", " << qubits_1 << "]"
		  << ", Quest = [" << quest_0 << ", " << quest_1 << "]"
		  << std::endl
		  << "Qubits = [" << distance_0 << ", " << distance_1 << "]"
		  << ", Qubits norms = [" << norm_0 << ", "
		  << norm_1 << "]" << std::endl;
    }
}


/// Implementation for one qubit gates with no arguments
template<typename F, typename G>
void impl(OneQubitGate<F> fn, OneQubitGate<G> gn, F & qubits, G & quest)
{
    unsigned nqubits = qubits.getNumQubits();
    
    // Loop through every pair of qubits
    for (unsigned n = 0; n < nqubits; n++) {
	std::cout << "Qubit no " << n << ": ";
	
	// Apply the gates stored in member pointers fn and gn
	std::invoke(fn, qubits, n);
	std::invoke(gn, quest, n);

	///\todo Should be float?
	double distance = fubiniStudy(qubits.getState(),
				      quest.getState());

	std::cout << "Distance = " << distance << std::endl;
    }
}

/// Implementation for single qubit gates with one double argument
template<typename F, typename G>
void impl(OneQubitGate<F, typename F::Fp_type> fn,
	  OneQubitGate<G, typename G::Fp_type> gn,
	  F & qubits, G & quest)
{
    unsigned nqubits = qubits.getNumQubits();
    std::vector<typename F::Fp_type> F_phase_list
	= makeRandomPhases<typename F::Fp_type>(nqubits);
    std::vector<typename G::Fp_type> G_phase_list
	= convertVector<typename G::Fp_type>(F_phase_list);
    
    // Loop through every pair of qubits
    for (unsigned n = 0; n < nqubits; n++) {
	std::cout << "Qubit no " << n << ": ";
	// Apply the gates stored in member pointers fn and gn
	std::invoke(fn, qubits, n, F_phase_list[n]);
	std::invoke(gn, quest, n, G_phase_list[n]);

	///\todo Should be float?
	double distance = fubiniStudy(qubits.getState(),
				      quest.getState());

	std::cout << "Distance = " << distance << std::endl;
    }
}

/// Implementation for two qubit gates with no arguments
template<typename F, typename G>
void impl(TwoQubitGate<F> fn, TwoQubitGate<G> gn,
	  F & qubits, G & quest)
{
    unsigned nqubits = qubits.getNumQubits();
    
    // Loop through every pair of qubits
    for (unsigned n = 0; n < nqubits; n++) {
	for (unsigned m = 0; m < nqubits; m++) {
	    if (n != m) {
		std::cout << "Qubits " << n  << " and " << m << ": ";

		// Apply the gates stored in member pointers fn and gn
		std::invoke(fn, qubits, n, m);
		std::invoke(gn, quest, n, m);
		
		double distance = fubiniStudy(qubits.getState(),
					      quest.getState());
		std::cout << "Distance = " << distance << std::endl;
	    }
	}
    }
}

/// Implementation for two qubit gates with one double argument
template<typename F, typename G>
void impl(TwoQubitGate<F, typename F::Fp_type> fn,
	  TwoQubitGate<G, typename G::Fp_type > gn,
	  F & qubits, G & quest)
{

    unsigned nqubits = qubits.getNumQubits();
    //std::vector<double> phase_list = makeRandomPhases(nqubits*nqubits);
    std::vector<typename F::Fp_type> F_phase_list
	= makeRandomPhases<typename F::Fp_type>(nqubits*nqubits);
    std::vector<typename G::Fp_type> G_phase_list
	= convertVector<typename G::Fp_type>(F_phase_list);
    
    // Loop through every pair of qubits
    for (unsigned n = 0; n < nqubits; n++) {
	for (unsigned m = 0; m < nqubits; m++) {
	    if (n != m) {
		std::cout << "Qubits " << n  << " and " << m << ": ";

		// Apply the gates stored in member pointers fn and gn
		std::invoke(fn, qubits, n, m, F_phase_list[nqubits*n + m]);
		std::invoke(gn, quest, n, m, G_phase_list[nqubits*n + m]);

		///\todo Should be float or Fp_type?
		double distance = fubiniStudy(qubits.getState(),
					      quest.getState());
		std::cout << "Distance = " << distance << std::endl;
	    }
	}
    }
}

template<typename F, typename G, typename... ArgsF, typename... ArgsG>
void verifyGate(Gate<F, ArgsF...> fn, Gate<G, ArgsG...> gn,
		unsigned nqubits, std::string)
{
    // Check for equal length argument lists
    ///\todo Need to do stronger check for compatible lists
    static_assert(sizeof...(ArgsF) == sizeof...(ArgsG),
		  "Cannot compare two gates with different argument lists");

    
    std::cout << "============================================" << std::endl;
    std::cout << "Verifying gates on " << nqubits << " qubits."
	      << std::endl;

    // Make a random state
    std::vector<complex<typename F::Fp_type>> state
	= makeRandomState<typename F::Fp_type>(nqubits);
    F qubits(state);
    G quest(convertState<typename G::Fp_type>(state));
	
    // Run the verification
    impl(fn, gn, qubits, quest);
}
