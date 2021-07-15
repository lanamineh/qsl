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

#ifndef QUBITS_NP_HPP
#define QUBITS_NP_HPP

/**
 * \defgroup qubits_constructors_np Constructors
 * \brief Constructors for the default simulator object
 */

/**
 * \defgroup qubits_gates_np Quantum gates
 * \brief One- and two- qubit gates in the Qubits class
 */

/**
 * \defgroup qubits_meas_np Measurement
 * \brief Functions related to measurement, probabilities, etc.
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

    public:

	const static std::string name; ///< The name of this simulator
    
	/// Expose Fp to users of the class
	using Fp_type = Fp; ///< \todo Is there a less bad way to do this?

	/**
	 * \brief Initialise the class with a fixed number of qubits
	 * \ingroup qubits_constructors_np
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
	 * \brief  Initialise the class from a pre-prepared state vector
	 * \ingroup qubits_constructors_np
	 *
	 * This function constructs an object with the given initial
	 * state vector. The state vector must have a length which 
	 * is a power of two, otherwise the function will throw 
	 * std::logic_error. Function also checks whether the state is
	 * number preserving and sets the corresponding number of ones.
	 *
	 * \param state A vector containing the initial state for the object
	 *
	 */
	Qubits(const std::vector<complex<Fp>> & state);

	/**
	 * \brief Copy constructor
	 * \ingroup qubits_constructors_np
	 *
	 * You can make copies of this object by constructing from an 
	 * object of the same type.
	 */
	Qubits(const Qubits & ) = default;

	/**
	 * \brief Copy-assignment operator
	 * \ingroup qubits_constructors_np
	 *
	 * You can assign one Qubits object to another, provided that they both
	 * represent the same number of qubits. In this case, this operation
	 * copies the state vector of one object to the other. If the number of
	 * qubits are not equal, this function throws std::logic_error.
	 */
	void operator = (const Qubits & old);

	/// Return the number of qubits
	unsigned getNumQubits() const;

	/// Return the number of ones
	unsigned getNumOnes() const;

	/// Set the number of ones, resets the state to lowest computational basis state.
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

	/// Initialise the lookup tables
	void initLookup();
    
	/**
	 * \brief Apply a phase shift to qubit number targ.
	 *
	 * \ingroup qubits_gates_np
	 *
	 * \f[ 
	 * R_\theta = \begin{pmatrix}
	 *            1 & 0 \\
	 *            0 & e^{i\theta} \\
	 *            \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to phase shift the qubit by.
	 */
	void phase(unsigned targ, Fp angle);

	/**
	 * \brief Apply the Pauli Z gate to qubit number targ.
	 *
	 * \ingroup qubits_gates_np
	 *
	 * \f[ 
	 * Z = \begin{pmatrix}
	 *     1 & 0 \\
	 *     0 & -1 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 */
	void pauliZ(unsigned targ);

	
	/**
	 * \brief Rotate around the z-axis of the Bloch sphere \f$ e^{-i\theta Z/2} \f$
	 *
	 * \ingroup qubits_gates_np
	 *
	 * This single qubit gate applies the following 2x2 matrix to each
	 * pair of \f$ |0\rangle \f$ and \f$ |1\rangle \f$ amplitudes for 
	 * angle \f$ \theta \f$:
	 * 
	 * \f[ 
	 * R_z = \begin{pmatrix}
	 *       e^{-i\theta/2} & 0 \\
	 *       0 & e^{i\theta/2} \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by.
	 */
	void rotateZ(unsigned targ, Fp angle);
	
	/**
	 * \brief Perform a controlled phase shift on two qubits. 
	 *
	 * \ingroup qubits_gates_np
	 *
	 * A phase is added if both qubits are in the \f$ |1\rangle \f$ state. 
	 *
	 * \f[ 
	 * CR_\theta = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 0 & 0 & e^{i\theta}
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, phase shift is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to phase shift the target qubit by.
	 */
	void controlPhase(unsigned ctrl, unsigned targ, Fp angle);


	/**
	 * \brief Perform a swap gate on two qubits. 
	 *
	 * \ingroup qubits_gates_np

	 * \f[ 
	 * SWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param q1 The first qubit to swap.
	 * \param q2 The second qubit to swap.
	 */
	void swap(unsigned q1, unsigned q2);

	/**
	 * \brief Perform a controlled Z gate on two qubits. 
	 *
	 * \ingroup qubits_gates_np
	 *
	 *
	 * \f[ 
	 * CZ = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 0 & 0 & -1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, Z is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void controlZ(unsigned ctrl, unsigned targ);

    
	/**
	 * \brief Measure a qubit and collapse the state to its outcome.
	 *
	 * \ingroup qubits_meas_np
	 *
	 * This is not a reversible operation unlike applying quantum gates.
	 *
	 * A random number is generated to determine whether the given qubit
	 * targ is measured to be 0 or 1. The state vector is then collapsed
	 * to that outcome by zeroing out all the amplitudes that do not
	 * correspond to the generated outcome. Note that the state vector
	 * does not change size.
	 *
	 * \param targ The qubit to measure.
	 * \return The value of the measured qubit (0 or 1).
	 */
	int measure(unsigned targ);

	/**
	 * \brief Measure all of the qubits at once and collapse to
	 * the resulting computational basis state.
	 *
	 * \ingroup qubits_meas_np
	 *
	 * Measuring all of the qubits at once is the same as measuring them 
	 * one by one. 
	 *
	 * \return The result of the measurement.
	 */
	std::size_t measureAll();
    
	/**
	 * \brief Calculate the probability of qubit targ being measured 
	 * in the given outcome (0 or 1).
	 *
	 * \ingroup qubits_meas_np
	 *
	 * \param targ The qubit to calculate the probability for.
	 * \param outcome The outcome (0 or 1) we are calculating the probability of.
	 * \return The probability of the qubit being measured in the given outcome.
	 */
	Fp prob(unsigned targ, unsigned outcome) const;

	/**
	 * \brief Perform a post-selection measurement. The state is collapsed 
	 * to the given outcome for the given qubit. 
	 *
	 * \ingroup qubits_meas_np
	 *
	 * This is not a reversible operation unlike applying quantum gates.
	 * The state vector is collapsed by zeroing out all the amplitudes 
	 * that do not correspond to the given outcome. Note that the state vector
	 * does not change size.
	 * 
	 * \param targ The qubit to measure.
	 * \param outcome The outcome (0 or 1) to post select on.
	 * \return The probability of measuring qubit targ in the given outcome.
	 */
	Fp postselect(unsigned targ, unsigned outcome);

	/**
	 * \brief Sample measurement outcome for one qubit multiple times.
	 * 
	 * \ingroup qubits_meas_np
	 *
	 * \param targ The qubit to measure.
	 * \param nsamples The number of samples to draw.
	 * \return A vector containing all the measurement results.
	 */
	std::vector<std::size_t> sample(unsigned targ, std::size_t nsamples);
    
	/**
	 * \brief Sample measurement outcome for all of the qubits 
	 * at once multiple times.
	 *
	 * \ingroup qubits_meas_np
	 *
	 * Measuring all of the qubits at once is the same as measuring them 
	 * one by one. This function implements a very efficient way of
	 * measuring all the qubits multiple times and should be used 
	 * instead of the measure function where possible. Note that this 
	 * function does not modify the state vector.
	 * 
	 * \param nsamples The number of measurements to perform.
	 * \return A map containing all of the measurement results, 
	 *         mapping outcomes to the number of times they happen. 
	 */
	std::map<std::size_t, std::size_t> sampleAll(std::size_t nsamples);

	/**
	 * \brief Sample measurement outcome for all of the qubits 
	 * at once multiple times.
	 *
	 * \ingroup qubits_meas_np
	 *
	 * Testing implementing sampling using the method in qsim. 
	 *
	 * \param nmeas The number of measurements to perform.
	 * \return A map containing all of the measurement results, 
	 *         mapping outcomes to the number of times they happen. 
	 */
	std::map<std::size_t, std::size_t> sampleAll2(std::size_t nsamples);
    };

}
    
#endif
