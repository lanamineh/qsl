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
 * \brief Checkers for use with the Verify class
 *
 */

#ifndef CHECKERS_HPP
#define CHECKERS_HPP

#include <memory>

/**
 * \brief Check gates
 *
 * This class contains the routines that are used to check 
 * gates. There is no checkAll method, so it cannot be used
 * as an argument in the Verify class. However, classes 
 * derived from this one that implement checkAll can be used
 * in the Verify class.
 *
 * The class does nothing other than implement gates. It does
 * not modify the simulators in any other way (e.g. setting 
 * states, etc.)
 *
 * \todo Find a way to control the phases, etc.
 */
template<typename Sim1, typename Sim2>
class GateChecker
{
    std::unique_ptr<Sim1> sim1;
    std::unique_ptr<Sim2> sim2;

public:
    
    /// Implementation for one qubit gates with no arguments
    void check(OneQubitGate<Sim1> fn, OneQubitGate<Sim2> gn)
	{
	    unsigned nqubits = sim1->getNumQubits();
	    
	    // Loop through every pair of sim1
	    for (unsigned n = 0; n < nqubits; n++) {
		std::cout << "Qubit no " << n << ": ";
	    
		// Apply the gates stored in member pointers fn and gn
		std::invoke(fn, sim1, n);
		std::invoke(gn, sim2, n);
	    
		///\todo Should be float?
		double distance = fubiniStudy(sim1->getState(),
					      sim2->getState());

		std::cout << "Distance = " << distance << std::endl;
	    }
	}

/// Checkementation for single qubit gates with one double argument
    void check(OneQubitGate<Sim1, typename Sim1::Fp_type> fn,
	       OneQubitGate<Sim2, typename Sim2::Fp_type> gn)
	{
	    unsigned nqubits = sim1->getNumQubits();
	    std::vector<typename Sim1::Fp_type> sim1_phase_list
		= makeRandomPhases<typename Sim1::Fp_type>(nqubits);
	    std::vector<typename Sim2::Fp_type> sim2_phase_list
		= convertVector<typename Sim2::Fp_type>(sim1_phase_list);
    
	    // Loop through every pair of sim1
	    for (unsigned n = 0; n < nqubits; n++) {
		std::cout << "Qubit no " << n << ": ";
		// Apply the gates stored in member pointers fn and gn
		std::invoke(fn, sim1, n, sim1_phase_list[n]);
		std::invoke(gn, sim2, n, sim2_phase_list[n]);

		///\todo Should be float?
		double distance = fubiniStudy(sim1->getState(),
					      sim2->getState());
		std::cout << "Distance = " << distance << std::endl;
	    }
	}

/// Checkementation for two qubit gates with no arguments
    void check(TwoQubitGate<Sim1> fn, TwoQubitGate<Sim2> gn)
	{
	    unsigned nqubits = sim1->getNumQubits();
    
	    // Loop through every pair of sim1
	    for (unsigned n = 0; n < nqubits; n++) {
		for (unsigned m = 0; m < nqubits; m++) {
		    if (n != m) {
			std::cout << "Sim1 " << n  << " and " << m << ": ";

			// Apply the gates stored in member pointers fn and gn
			std::invoke(fn, sim1, n, m);
			std::invoke(gn, sim2, n, m);
		
			double distance = fubiniStudy(sim1->getState(),
						      sim2->getState());
			std::cout << "Distance = " << distance << std::endl;
		    }
		}
	    }
	}

/// Checkementation for two qubit gates with one double argument
    void check(TwoQubitGate<Sim1, typename Sim1::Fp_type> fn,
	       TwoQubitGate<Sim2, typename Sim2::Fp_type> gn)
	{
	
	    unsigned nqubits = sim1->getNumQubits();
	    //std::vector<double> phase_list = makeRandomPhases(nqubits*nqubits);
	    std::vector<typename Sim1::Fp_type> sim1_phase_list
		= makeRandomPhases<typename Sim1::Fp_type>(nqubits*nqubits);
	    std::vector<typename Sim2::Fp_type> sim2_phase_list
		= convertVector<typename Sim2::Fp_type>(sim1_phase_list);
    
	    // Loop through every pair of sim1
	    for (unsigned n = 0; n < nqubits; n++) {
		for (unsigned m = 0; m < nqubits; m++) {
		    if (n != m) {
			std::cout << "Sim1 " << n  << " and " << m << ": ";

			// Apply the gates stored in member pointers fn and gn
			std::invoke(fn, sim1, n, m, sim1_phase_list[nqubits*n + m]);
			std::invoke(gn, sim2, n, m, sim2_phase_list[nqubits*n + m]);

			///\todo Should be float or Fp_type?
			double distance = fubiniStudy(sim1->getState(),
						      sim2->getState());
			std::cout << "Distance = " << distance << std::endl;
		    }
		}
	    }
	}


    /**
     * \brief Attach simulator objects to this checker
     */
    void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	{
	    std::swap(sim1, p1);
	    std::swap(sim2, p2);
	}

    void configureChecker(int val)
	{
	    std::cout << "Doing something different: " << val << std::endl;
	}
    
}; 

/// Example gate checker (check all gates)
template<HasAllGates Sim1, HasAllGates Sim2>
class DefaultGateChecker : public GateChecker<Sim1,Sim2>
{
public:
    void checkAll() {
	std::cout << "Checking pauliX" << std::endl;
	this->check(&Sim1::pauliX, &Sim2::pauliX);
	std::cout << "Checking rotateX" << std::endl;
	this->check(&Sim1::rotateX, &Sim2::rotateX);
	std::cout << "Checking phase" << std::endl;
	this->check(&Sim1::phase, &Sim2::phase);
	std::cout << "Checking controlNot" << std::endl;
	this->check(&Sim1::controlNot, &Sim2::controlNot);
	std::cout << "Checking controlPhase" << std::endl;
	this->check(&Sim1::controlPhase, &Sim2::controlPhase);
    }
};

/**
 * \brief Check number preserved gates
 *
 * Pass this class to the Verify class to check number
 * preserving gates.
 *
 */
template<HasNPGates Sim1, HasNPGates Sim2>
class NPGateChecker : public GateChecker<Sim1,Sim2>
{
public:
    void checkAll() {
	std::cout << "Checking phase" << std::endl;
	this->check(&Sim1::phase, &Sim2::phase);
	std::cout << "Checking controlPhase" << std::endl;
	this->check(&Sim1::controlPhase, &Sim2::controlPhase);
    }
};

/**
 * \brief Check measurement related functions
 *
 * The purpose of this class is to contain checks for all the
 * measurement functions. It is not supposed to be used directly
 * (it doesn't have a checkAll function). Instead, classes
 * derived from this class can call a subset of the measurement 
 * checks depending on what is required.
 *
 */
template<typename Sim1, typename Sim2>
class PostselectChecker
{
    std::unique_ptr<Sim1> sim1;
    std::unique_ptr<Sim2> sim2;

public:

    /**
     * \brief Check postselection
     *
     * For concepts: requires copy assignment, copy, postselect and getState.
     */
    void checkAll()
	{
	    if (not sim1 or not sim2) {
		throw std::logic_error(
		    "Cannot run check before binding simulators");
	    }
	    
	    unsigned nqubits = sim1->getNumQubits();
	    std::cout << "============================================"
		      << std::endl;
	    std::cout << "Verifying postselect on " << nqubits << " qubits."
		      << std::endl;

	    std::cout << "Printing probabilities of collapsing to outcome"
		      << std::endl;
    
	    for (unsigned n = 0; n < nqubits; n++) {

		// Make local copy
		Sim1 q1(*sim1); // Make a local copy
		Sim2 q2(*sim2); // Make a local copy
	
		// Postselect on 0
		double qubits_0 = q1.postselect(n, 0);
		double quest_0 = q2.postselect(n, 0);
		double norm_0 = norm(q1.getState());
		double distance_0 = fubiniStudy(q1.getState(), q2.getState());

		// Reset the state vectors
		q1 = *sim1;
		q2 = *sim2;

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

    /**
     * \brief Attach simulator objects to this checker
     */
    void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	{
	    std::swap(sim1, p1);
	    std::swap(sim2, p2);
	}
    
    
};

#endif
