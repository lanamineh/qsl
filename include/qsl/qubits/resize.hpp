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

#ifndef QUBITS_RESIZE_HPP
#define QUBITS_RESIZE_HPP

/**
 * \defgroup qubits_constructors Constructors
 * \brief Constructors for the resizeable simulator object
 */

/**
 * \defgroup qubits_gates Quantum gates
 * \brief One- and two- qubit gates in the Qubits class
 */

/**
 * \defgroup qubits_meas Measurement
 * \brief Functions related to measurement, probabilities, etc.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <map>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"

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
	 * \ingroup qubits_constructors
	 *
	 * This function constructors an object with the specified 
	 * number of qubits. The simulator is initialised in the
	 * all zero state.
	 *
	 * \param nqubits The number of qubits to simulate
	 */ 
	Qubits(unsigned nqubits);

	/**
	 * \brief  Initialise the class from a pre-prepared state vector
	 * \ingroup qubits_constructors
	 *
	 * This function constructs an object with the given initial
	 * state vector. The  state vector must have a length which 
	 * is a power of two, otherwise the function will throw 
	 * std::logic_error.
	 *
	 * \param state A vector containing the initial state for the object
	 *
	 */
	Qubits(const std::vector<complex<Fp>> & state);

	/**
	 * \brief Copy constructor
	 * \ingroup qubits_constructors
	 *
	 * You can make copies of this object by constructing from an 
	 * object of the same type.
	 */
	Qubits(const Qubits & ) = default;

	/**
	 * \brief Copy-assignment operator
	 * \ingroup qubits_constructors
	 *
	 * You can assign one Qubits object to another, provided that they both
	 * represent the same number of qubits. In this case, this operation
	 * copies the state vector of one object to the other. If the number of
	 * qubits are not equal, this function throws std::logic_error.
	 */
	void operator = (const Qubits & old);

	/// Return the number of qubits
	unsigned getNumQubits() const;
    
	/// Reset to the all-zero state
	void reset();

	/// Print the state vector
	void print() const;

	/// Get the state vector associated to the qubits
	std::vector<complex<Fp>> getState() const;

	/// Set the state vector (i.e. re-initialise the state vector)
	void setState(const std::vector<complex<Fp>> & state);

	/// Set the state vector to a computational basis state
	void setBasisState(std::size_t index);

	/**
	 * \brief Add a qubit to the state vector
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
	 * \brief Rotate around the x-axis of the Bloch sphere \f$ e^{-i\theta X/2} \f$
	 *
	 * \ingroup qubits_gates
	 *
	 * This single qubit gate applies the following 2x2 matrix to each
	 * pair of \f$ |0\rangle \f$ and \f$ |1\rangle \f$ amplitudes for 
	 * angle \f$ \theta \f$:
	 * 
	 * \f[ 
	 * R_x = \begin{pmatrix}
	 *       \cos(\theta/2) & -i\sin(\theta/2) \\
	 *       -i\sin(\theta/2) & \cos(\theta/2) \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by.
	 */
	void rotateX(unsigned targ, Fp angle);

	/**
	 * \brief Apply the Hadamard gate to qubit number targ.
	 *
	 * \ingroup qubits_gates
	 *
	 * \f[ 
	 * H = \frac{1}{\sqrt{2}}\begin{pmatrix}
	 *     1 & 1 \\
	 *     1 & -1 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 */
	void hadamard(unsigned targ);
	
	/**
	 * \brief Apply the Pauli X gate to qubit number targ.
	 *
	 * \ingroup qubits_gates
	 *
	 * \f[ 
	 * X = \begin{pmatrix}
	 *     0 & 1 \\
	 *     1 & 0 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 */
	void pauliX(unsigned targ);

	/**
	 * \brief Apply a phase shift to qubit number targ.
	 *
	 * \ingroup qubits_gates
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
	 * \brief Perform the controlled Not (CNOT) gate on two qubits.
	 *
	 * \ingroup qubits_gates
	 *
	 * Controlled on the first qubit, the matrix is:
	 * \f[ 
	 * CNOT = \begin{pmatrix}
	 *        1 & 0 & 0 & 0 \\
	 *        0 & 1 & 0 & 0 \\
	 *        0 & 0 & 0 & 1 \\
	 *        0 & 0 & 1 & 0
	 *        \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, NOT is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void controlNot(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled phase shift on two qubits. 
	 *
	 * \ingroup qubits_gates
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
	 * \ingroup qubits_gates

	 * \f[ 
	 * CR_\theta = \begin{pmatrix}
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
	 * \brief Measure a qubit and collapse the state to its outcome.
	 *
	 * \ingroup qubits_meas
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
	 * \brief Measure a qubit and remove the qubit from the state vector.
	 *
	 * \ingroup qubits_meas
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
	 *
	 * \ingroup qubits_meas
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
	 * \ingroup qubits_meas
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
	 * \ingroup qubits_meas
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
	 * \ingroup qubits_meas
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
	 * \ingroup qubits_meas
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
	 * \ingroup qubits_meas
	 *
	 * Testing implementing sampling using the method in qsim. 
	 *
	 * \param nmeas The number of measurements to perform.
	 * \return A map containing all of the measurement results, 
	 *         mapping outcomes to the number of times they happen. 
	 */
	std::map<std::size_t, std::size_t> sampleAll2(std::size_t nsamples);
    };

    
    // Explicit instantiation declarations are required to avoid
    // compiler warnings in clang, when template instantiations
    // appear in another translation unit
#ifdef __clang__
    template<> const std::string Qubits<Type::Default, double>::name;
    template<> const std::string Qubits<Type::Default, float>::name;
#endif

}
    
#endif