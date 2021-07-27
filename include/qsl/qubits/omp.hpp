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
 * \file qubits/omp.hpp
 * \brief Contains the OpenMP simulator object for manipulating qubits
 * 
 */

#ifndef QSL_QUBITS_OMP_HPP
#define QSL_QUBITS_OMP_HPP

/**
 * \defgroup qubits_omp functionsOmp
 * \brief Functions that are different/specific to the OMP simulator.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"
#include <omp.h>
#include <iostream>

namespace qsl {

    /**
     * \brief OpenMP quantum simulator object
     *
     * This classes implements parallelism using OpenMP. Otherwise, the
     * implementation is very similar to Qubits<Type::Default, Fp>. You
     * can specify the number of threads to use as an optional parameter
     * in the constructor.
     *
     */
    template<std::floating_point Fp>
    class Qubits<Type::Omp, Fp>
    {
	const unsigned nqubits;
	const std::size_t dim; ///< The length of the state vector
	const unsigned nthreads; ///< The number of threads
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
	 * \brief Generate the cumulative probability vector, ignoring
	 * amplitudes that are zero.
	 */
	std::vector<Dist> generateDist();

    
	void pauliX_inner(std::size_t start, std::size_t end, std::size_t k);
	void pauliX_outer(std::size_t start, std::size_t end, std::size_t k);
    
    public:

	const static std::string name; ///< The name of this simulator
    
	/// Expose Fp to users of the class
	using Fp_type = Fp; ///< \todo Is there a less bad way to do this?
    
	/**
	 * \brief Initialise the class with a fixed number of qubits
	 * \ingroup qubits_omp
	 *
	 * This function constructors an object with the specified 
	 * number of qubits and number of threads. The simulator
	 * is initialised in the all zero state.
	 *
	 * \param nqubits The number of qubits to simulate
	 * \param nthreads The number of threads to use
	 */ 
	Qubits(unsigned nqubits, unsigned nthreads = 4);

	/**
	 * \brief  Initialise the class from a pre-prepared state vector
	 * \ingroup qubits_omp
	 *
	 * This function constructs an object with the given initial
	 * state vector, and the specified number of threads. The 
	 * state vector must have a length which is a power of two,
	 * otherwise the function will throw std::logic_error.
	 *
	 * \param state A vector containing the initial state for the object
	 * \param nthreads The number of threads to use
	 */
	Qubits(const std::vector<complex<Fp>> & state, unsigned nthreads = 4);

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
	 * \brief Sample measurement outcome for one qubit multiple times.
	 */
	std::vector<std::size_t> sample(unsigned targ, std::size_t nsamples);
    
	/**
	 * \brief Sample measurement outcome for all of the qubits 
	 * at once multiple times.
	 */
	std::map<std::size_t, std::size_t> sampleAll(std::size_t nsamples);

    };

}
    
#endif
