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

#ifndef QUBITS_OMP_HPP
#define QUBITS_OMP_HPP

/**
 * \defgroup qubits_omp_constructors Constructors
 * \brief Constructors for the OpenMP simulator object
 */

/**
 * \defgroup qubits_omp_gates Quantum gates
 * \brief One- and two- qubit gates in the Qubits class
 */

/**
 * \defgroup qubits_meas Measurement
 * \brief Functions related to measurement, probabilities, etc.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"
#include <omp.h>

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
	 * \ingroup qubits_omp_constructors
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
	 * \ingroup qubits_omp_constructors
	 *
	 * This function constructs an object with the given initial
	 * state vector, and the specified number of threads. The 
	 * state vector must have a length which is a power of two,
	 * otherwise the function will throw std::logic_error.
	 *
	 * \param state A vector containing the initial state for the object
	 * \param nthreads The number of threads to use
	 *
	 */
	Qubits(const std::vector<complex<Fp>> & state, unsigned nthreads = 4);

	/**
	 * \brief Copy constructor
	 * \ingroup qubits_omp_constructors
	 *
	 * You can make copies of this object by constructing from an 
	 * object of the same type. 
	 */
	Qubits(const Qubits & ) = default;
    
	/**
	 * \brief Copy-assignment operator
	 * \ingroup qubits_omp_constructors
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
	 * \brief Rotate around the x-axis of the Bloch sphere 
	 * \f$ e^{-i\theta X/2} \f$
	 * \ingroup qubits_omp_gates
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
	 * \brief Apply the Pauli X gate to qubit number targ.
	 * \ingroup qubits_omp_gates
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
	 * \ingroup qubits_omp_gates
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
	 * \ingroup qubits_omp_gates
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
	 * \ingroup qubits_omp_gates
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
	 * \brief Measure a qubit and collapse the state to its outcome.
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
	 * \brief Measure all of the qubits at once and collapse to
	 * the resulting computational basis state.
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
	double prob(unsigned targ, unsigned outcome) const;

	/**
	 * \brief Perform a post-selection measurement. The state is collapsed 
	 * to the given outcome for the given qubit. 
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
	double postselect(unsigned targ, unsigned outcome);

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

    };

}
    
#endif
