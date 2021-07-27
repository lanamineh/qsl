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
 * \file qubits/np.hpp
 * \brief Contains the number preserving simulator.
 * 
 */

#ifndef QSL_QUBITS_NP_HPP
#define QSL_QUBITS_NP_HPP

/**
 * \defgroup qubits_np functionsNP
 * \brief Functions that are different/specific to the NP simulator.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <map>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"
#include <iostream>

namespace qsl {

    /**
     * \brief Numper preserving simulator object
     *
     * This version of the Qubits object only has non-zero amplitudes
     * associated to bit strings with a specified number of ones.
     *
     */
    template<std::floating_point Fp>
    class Qubits<Type::NP, Fp>
    {
	const unsigned nqubits;
	const std::size_t dim; ///< The length of the state vector
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
	 * \ingroup qubits_np
	 *
	 * This function constructors an object with the specified 
	 * number of qubits. The simulator is initialised in the
	 * all zero state.
	 *
	 * \param nqubits The number of qubits to simulate
	 * \param nones The number of ones to preserved throughout
	 */ 
	Qubits(unsigned nqubits, unsigned nones = 1);

	/**
	 * \brief Initialise the class from a pre-prepared state vector
	 *
	 * Function also checks whether the state is
	 * number preserving and sets the corresponding number of ones.
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

	/// Return the number of ones
	unsigned getNumOnes() const;

	/**
	 * \brief Set the number of ones, resets the state to lowest 
	 * computational basis state.
	 * \ingroup qubits_np
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
