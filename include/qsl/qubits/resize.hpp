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
 * \file qubits/resize.hpp
 * \brief Contains the resizeable variant of the Qubits object
 * 
 */

/**
 * \defgroup qubits_resize functionsResize
 * \brief Functions that are different/specific to the Resize simulator.
 */


#ifndef QSL_QUBITS_RESIZE_HPP
#define QSL_QUBITS_RESIZE_HPP

#include <cstddef>
#include <cstdint>
#include <vector>
#include <map>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"
#include <iostream>

namespace qsl {
    
    /**
     * \brief Default quantum simulator object
     *
     * This version of the Qubits class does not use threading
     * or any inline assembly of any kind.
     *
     */
    template<std::floating_point Fp>
    class Qubits<Type::Resize, Fp>
    {
	unsigned nqubits;
	std::size_t dim; ///< The length of the state vector
	std::vector<complex<Fp>> state; ///< State vector for the qubits

	qsl::Random<Fp> random;

	struct Dist {
	    std::size_t index;
	    double prob;
	};
    
	/**
	 * \brief Run binary search to draw a sample from the cumulative
	 * probability vector.
	 */
	std::size_t drawSample(const std::vector<Dist> & dist);

	/**
	 * \brief Collapse a qubit to the given outcome with the given
	 * renormalisation factor.
	 */
	void collapse(unsigned targ, unsigned outcome, Fp factor);

	/**
	 * \brief Collapse a qubit to the given outcome and remove
	 * from the state vector.
	 * \ingroup qubits_resize
	 * 
	 * This member function has the same behaviour as the collapse
	 * function, but the measured qubit is also removed from the state
	 * vector, and the number of qubits returned by getNumQubits() is reduced
	 * by one. 
	 *
	 * The resulting state vector preserves the order of the other qubits,
	 * meaning that a qubit at location index > targ will now be addressable
	 * using index - 1. 
	 *
	 */
	void collapseOut(unsigned targ, unsigned outcome, Fp factor);
	
	/**
	 * \brief Generate the cumulative probability vector, ignoring
	 * amplitudes that are zero.
	 */
	std::vector<Dist> generateDist();
    
    public:

	const static std::string name; ///< The name of this simulator
    
	/// Expose Fp to users of the class
	using Fp_type = Fp; ///< \todo Is there a less bad way to do this?

	/**
	 * \brief Initialise the class with a fixed number of qubits
	 */ 
	Qubits(unsigned nqubits);

	/**
	 * \brief  Initialise the class from a pre-prepared state vector
	 */
	Qubits(const std::vector<complex<Fp>> & state);

	/**
	 * \brief Copy constructor
	 */
	Qubits(const Qubits & ) = default;

	/**
	 * \brief Copy-assignment operator
	 */
	void operator = (const Qubits & old);

	/// Return the number of qubits
	unsigned getNumQubits() const;
    
	/// Reset to the all-zero state
	void reset();

	/// Print the state vector
	void print(std::ostream & os = std::cout) const;

	/// Get the state vector associated to the qubits
	std::vector<complex<Fp>> getState() const;

	/// Set the state vector (i.e. re-initialise the state vector)
	void setState(const std::vector<complex<Fp>> & state);

	/// Set the state vector to a computational basis state
	void setBasisState(std::size_t index);

	/**
	 * \brief Add a qubit to the state vector
	 * \ingroup qubits_resize
	 *
	 * This function adds a qubit to the end of the list of qubits.
	 * The resulting state is a tensor product between the previous
	 * state vector and the |0) state of the newly added qubit.
	 *
	 * The index of the new qubit is one larger than the index of
	 * the largest qubit in the previous state vector.
	 *
	 * Calling getNumQubits() will return the new number of qubits,
	 * which is one larger than before the function is called.
	 *
	 */
	void addQubit();

	/**
	 * \brief Rotate around the x-axis of the Bloch sphere 
	 */
	void rotateX(unsigned targ, Fp angle);

	/**
	 * \brief Rotate around the y-axis of the Bloch sphere 
	 */
	void rotateY(unsigned targ, Fp angle);
	
	/**
	 * \brief Rotate around the z-axis of the Bloch sphere
	 */
	void rotateZ(unsigned targ, Fp angle);
	
	/**
	 * \brief Apply the Hadamard gate to qubit number targ.
	 */
	void hadamard(unsigned targ);
	
	/**
	 * \brief Apply the Pauli X gate to qubit number targ.
	 */
	void pauliX(unsigned targ);

	/**
	 * \brief Apply the Pauli Y gate to qubit number targ.
	 */
	void pauliY(unsigned targ);

	/**
	 * \brief Apply the Pauli Z gate to qubit number targ.
	 */
	void pauliZ(unsigned targ);
	
	/**
	 * \brief Apply a phase shift to qubit number targ.
	 */
	void phase(unsigned targ, Fp angle);

	/**
	 * \brief Perform the controlled Not (CNOT) gate on two qubits.
	 */
	void controlNot(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform the CY gate on two qubits.
	 */
	void controlY(unsigned ctrl, unsigned targ);
	
	/**
	 * \brief Perform the CZ gate on two qubits.
	 */
	void controlZ(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform the CRx gate on two qubits.
	 */
	void controlRotateX(unsigned ctrl, unsigned targ, Fp angle);

	/**
	 * \brief Perform the CRy gate on two qubits.
	 */
	void controlRotateY(unsigned ctrl, unsigned targ, Fp angle);

	/**
	 * \brief Perform the CRz gate on two qubits.
	 */
	void controlRotateZ(unsigned ctrl, unsigned targ, Fp angle);
	
	/**
	 * \brief Perform a controlled phase shift on two qubits.
	 */
	void controlPhase(unsigned ctrl, unsigned targ, Fp angle);

	/**
	 * \brief Perform the CH gate on two qubits.
	 */
	void controlHadamard(unsigned ctrl, unsigned targ);
	
	/**
	 * \brief Perform a swap gate on two qubits. 
	 */
	void swap(unsigned q1, unsigned q2);

	/**
	 * \brief Perform a fermionic swap gate on two qubits. 
	 */
	void fswap(unsigned q1, unsigned q2);
	
	/**
	 * \brief Measure a qubit and collapse the state to its outcome.
	 */
	int measure(unsigned targ);

	/**
	 * \brief Measure a qubit and remove the qubit from the state vector.
	 * \ingroup qubits_resize
	 *
	 * This is not a reversible operation unlike applying quantum gates.
	 * 
	 * This member function has the same behaviour as the measure(unsigned)
	 * function, but the measured qubit is also removed from the state
	 * vector, and the number of qubits returned by getNumQubits() is reduced
	 * by one. 
	 *
	 * The resulting state vector preserves the order of the other qubits,
	 * meaning that a qubit at location index > targ will now be addressable
	 * using index - 1.
	 * 
	 * \param targ The qubit to measure and remove.
	 * \return The value of the measured qubit (0 or 1).
	 */
	int measureOut(unsigned targ);
	
	/**
	 * \brief Measure all of the qubits at once and collapse to
	 * the resulting computational basis state.
	 */
	std::size_t measureAll();
    
	/**
	 * \brief Calculate the probability of qubit targ being measured 
	 * in the given outcome (0 or 1).
	 */
	Fp prob(unsigned targ, unsigned outcome) const;

	/**
	 * \brief Perform a post-selection measurement. The state is collapsed 
	 * to the given outcome for the given qubit. 
	  */
	Fp postselect(unsigned targ, unsigned outcome);

	/**
	 * \brief Perform a postselection measurement and
	 * remove the qubit from the state vector.
	 * \ingroup qubits_resize	 
	 *
	 * This is not a reversible operation unlike applying quantum gates.
	 * 
	 * This member function has the same behaviour as the postselect
	 * function, but the measured qubit is also removed from the state
	 * vector, and the number of qubits returned by getNumQubits() is reduced
	 * by one. 
	 *
	 * The resulting state vector preserves the order of the other qubits,
	 * meaning that a qubit at location index > targ will now be addressable
	 * using index - 1.
	 *
	 * \param targ The qubit to measure.
	 * \param outcome The outcome (0 or 1) to post select on.
	 * \return The probability of measuring qubit targ in the given outcome.
	 */
	Fp postselectOut(unsigned targ, unsigned outcome);
	
	/**
	 * \brief Sample measurement outcome for one qubit multiple times.
	 */
	std::vector<std::size_t> sample(unsigned targ, std::size_t nsamples);
    
	/**
	 * \brief Sample measurement outcome for all of the qubits 
	 * at once multiple times.
	 */
	std::map<std::size_t, std::size_t> sampleAll(std::size_t nsamples);

	/**
	 * \brief Sample measurement outcome for all of the qubits 
	 * at once multiple times.
	 */
	std::map<std::size_t, std::size_t> sampleAll2(std::size_t nsamples);
    };
    
    // Explicit instantiation declarations are required to avoid
    // compiler warnings in clang, when template instantiations
    // appear in another translation unit
#ifdef __clang__
    template<> const std::string Qubits<Type::Resize, double>::name;
    template<> const std::string Qubits<Type::Resize, float>::name;
#endif

}
    
#endif
