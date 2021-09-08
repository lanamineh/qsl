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
 * \file qubits/omp_np.hpp
 * \brief Contains the number preserving simulator with multithreading.
 * 
 */

#ifndef QSL_QUBITS_OMP_NP_HPP
#define QSL_QUBITS_OMP_NP_HPP

/**
 * \defgroup qubits_omp_np functionsOmpNP
 * \brief Functions that are different/specific to the OmpNP simulator.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <map>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"
#include <iostream>
#include <omp.h>

namespace qsl {

    /**
     * \brief Numper preserving simulator object with OpenMP
     *
     * This version of the Qubits object only has non-zero amplitudes
     * associated to bit strings with a specified number of ones.
     *
     */
    template<std::floating_point Fp>
    class Qubits<Type::OmpNP, Fp>
    {
	const unsigned nqubits;
	const std::size_t dim; ///< The length of the state vector
	const unsigned nthreads; ///< The number of threads
	unsigned nones; ///< Number of ones to be preserved
	std::vector<complex<Fp>> state; ///< State vector for the qubits

	/// Set of lookup tables required for various gates 
	std::map<std::pair<unsigned, unsigned>, std::vector<std::size_t>> lookup;
    
	qsl::Random<Fp> random;

	struct Dist {
	    std::size_t index;
	    Fp prob;
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

	/// Initialise the lookup tables
	void initLookup();
	
    public:

	const static std::string name; ///< The name of this simulator
    
	/// Expose Fp to users of the class
	using Fp_type = Fp; ///< \todo Is there a less bad way to do this?

	/**
	 * \brief Initialise the class with a fixed number of qubits
	 * \ingroup qubits_omp_np
	 *
	 * This function constructors an object with the specified 
	 * number of qubits. The simulator is initialised in the
	 * all zero state.
	 *
	 * \param nqubits The number of qubits to simulate
	 * \param nones The number of ones to preserved throughout
	 */ 
	Qubits(unsigned nqubits, unsigned nones = 1, unsigned nthreads = 4);

	/**
	 * \brief Initialise the class from a pre-prepared state vector
	 * \ingroup qubits_omp_np
	 *
	 * Function also checks whether the state is
	 * number preserving and sets the corresponding number of ones.
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

	/// Return the number of ones
	unsigned getNumOnes() const;

	/**
	 * \brief Set the number of ones, resets the state to lowest 
	 * computational basis state.
	 */ 
	void setNumOnes(unsigned nones);
    
	/// Reset to the computational basis state with the lowest numerical value
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
	 * \brief Apply a phase shift to qubit number targ.
	 */
	void phase(unsigned targ, Fp angle);

	/**
	 * \brief Apply the Pauli Z gate to qubit number targ.
	 */
	void pauliZ(unsigned targ);

	/**
	 * \brief Rotate around the z-axis of the Bloch sphere 
	 */
	void rotateZ(unsigned targ, Fp angle);
	
	/**
	 * \brief Perform a controlled Z gate on two qubits. 
	 */
	void controlZ(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform the CRz gate on two qubits.
	 */
	void controlRotateZ(unsigned ctrl, unsigned targ, Fp angle);
	
	/**
	 * \brief Perform a controlled phase shift on two qubits. 
	 */
	void controlPhase(unsigned ctrl, unsigned targ, Fp angle);

	/**
	 * \brief Perform a swap gate on two qubits. 
	 */
	void swap(unsigned q1, unsigned q2);

	/**
	 * \brief Perform a fermionic swap gate on two qubits. 
	 */
	void fswap(unsigned q1, unsigned q2);

	/**
	 * \brief Perform a number preserved rotate X on two qubits. 
	 */
	void npRotateX(unsigned q1, unsigned q2, Fp angle);

	/**
	 * \brief Perform a number preserved rotate Y on two qubits. 
	 */
	void npRotateY(unsigned q1, unsigned q2, Fp angle);

	/**
	 * \brief Perform a number preserved Hadamard on two qubits. 
	 */
	void npHadamard(unsigned q1, unsigned q2);
	
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

	/**
	 * \brief Sample measurement outcome for all of the qubits 
	 * at once multiple times.
 	 */
	std::map<std::size_t, std::size_t> sampleAll2(std::size_t nsamples);
    };

}
    
#endif
