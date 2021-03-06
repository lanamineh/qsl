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
 * \file qubits/default.hpp
 * \brief Contains the default simulator object for manipulating qubits
 * 
 */

#ifndef QSL_QUBITS_DEFAULT_HPP
#define QSL_QUBITS_DEFAULT_HPP

/**
 * \defgroup qubits_constructors Constructors
 * \brief Constructors for the default simulator object
 */

/**
 * \defgroup qubits_utils Utilities
 * \brief Utilities for the default simulator object
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
    class Qubits<Type::Default, Fp>
    {
	const unsigned nqubits;
	const std::size_t dim; ///< The length of the state vector
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

	/**
	 * \brief Return the number of qubits
	 * \ingroup qubits_utils
	 */
	unsigned getNumQubits() const;
    
	/**
	 * \brief Reset to the all-zero state
	 * \ingroup qubits_utils
	 */
	void reset();

	/**
	 * \brief Print the state vector
	 * \ingroup qubits_utils
	 */
	void print(std::ostream & os = std::cout) const;

	/**
	 * \brief Get the state vector associated to the qubits
	 * \ingroup qubits_utils
	 * \todo Return std::vector<std::complex<Fp>> instead
	 */
	std::vector<complex<Fp>> getState() const;

	/**
	 * \brief Set the state vector (i.e. re-initialise the state vector)
	 * \ingroup qubits_utils
	 * \todo Take std::vector<std::complex<Fp>> instead
	 */
	void setState(const std::vector<complex<Fp>> & state);

	/**
	 * \brief Set the state vector to a computational basis state
	 * \ingroup qubits_utils
	 */
	void setBasisState(std::size_t index);
    
	/**
	 * \brief Rotate around the x-axis of the Bloch sphere 
	 * \f$ e^{-i\theta X/2} \f$
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
	 * \brief Rotate around the y-axis of the Bloch sphere 
	 * \f$ e^{-i\theta Y/2} \f$
	 *
	 * \ingroup qubits_gates
	 *
	 * This single qubit gate applies the following 2x2 matrix to each
	 * pair of \f$ |0\rangle \f$ and \f$ |1\rangle \f$ amplitudes for 
	 * angle \f$ \theta \f$:
	 * 
	 * \f[ 
	 * R_y = \begin{pmatrix}
	 *       \cos(\theta/2) & -\sin(\theta/2) \\
	 *       \sin(\theta/2) & \cos(\theta/2) \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by.
	 */
	void rotateY(unsigned targ, Fp angle);
	
	/**
	 * \brief Rotate around the z-axis of the Bloch sphere 
	 * \f$ e^{-i\theta Z/2} \f$
	 *
	 * \ingroup qubits_gates
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
	 * \brief Apply the Pauli Y gate to qubit number targ.
	 *
	 * \ingroup qubits_gates
	 *
	 * \f[ 
	 * Y = \begin{pmatrix}
	 *     0 & -i \\
	 *     i & 0 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * \param targ The target qubit.
	 */
	void pauliY(unsigned targ);
	
	/**
	 * \brief Apply the Pauli Z gate to qubit number targ.
	 *
	 * \ingroup qubits_gates
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
	 * \param ctrl The control qubit, X is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void controlNot(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled Y gate on two qubits. 
	 *
	 * \ingroup qubits_gates
	 *
	 *
	 * \f[ 
	 * CY = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & -i \\
	 *             0 & 0 & i & 0
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, Y is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void controlY(unsigned ctrl, unsigned targ);
	
	/**
	 * \brief Perform a controlled Z gate on two qubits. 
	 *
	 * \ingroup qubits_gates
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
	 * \brief Perform a controlled Rx gate on two qubits. 
	 *
	 * \ingroup qubits_gates
	 *
	 *
	 * \f[ 
	 * CRx = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & \cos(\theta/2) & -i\sin(\theta/2) \\
	 *             0 & 0 & -i\sin(\theta/2) & \cos(\theta/2)
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, Rx is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by.
	 */
	void controlRotateX(unsigned ctrl, unsigned targ, Fp angle);

	/**
	 * \brief Perform a controlled Ry gate on two qubits. 
	 *
	 * \ingroup qubits_gates
	 *
	 *
	 * \f[ 
	 * CRy = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & \cos(\theta/2) & -\sin(\theta/2) \\
	 *             0 & 0 & \sin(\theta/2) & \cos(\theta/2)
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, Ry is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by.
	 */
	void controlRotateY(unsigned ctrl, unsigned targ, Fp angle);
	
	/**
	 * \brief Perform a controlled Rz gate on two qubits. 
	 *
	 * \ingroup qubits_gates
	 *
	 *
	 * \f[ 
	 * CRz = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & e^{-i\theta/2} & 0 \\
	 *             0 & 0 & 0 & e^{i\theta/2}
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, Rz is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by.
	 */
	void controlRotateZ(unsigned ctrl, unsigned targ, Fp angle);
        
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
	 * \brief Perform a controlled H gate on two qubits. 
	 *
	 * \ingroup qubits_gates
	 *
	 *
	 * \f[ 
	 * CH = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 1/\sqrt{2} & 1/\sqrt{2} \\
	 *             0 & 0 & 1/\sqrt{2} & -1/\sqrt{2}
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param ctrl The control qubit, H is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void controlHadamard(unsigned ctrl, unsigned targ);
	
	/**
	 * \brief Perform a swap gate on two qubits. 
	 *
	 * \ingroup qubits_gates

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
	 * \brief Perform a fermionic swap gate on two qubits. 
	 *
	 * \ingroup qubits_gates

	 * \f[ 
	 * FSWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & -1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param q1 The first qubit to fswap.
	 * \param q2 The second qubit to fswap.
	 */
	void fswap(unsigned q1, unsigned q2);

	/**
	 * \brief Perform an X rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. This is equivalent to applying 
	 * \f$ e^{-i\theta (XX+YY)/2} \f$. 
	 *
	 * \ingroup qubits_gates

	 * \f[ 
	 *  npR_x(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -i\sin(\theta/2) & 0 \\
	 *             0 & -i\sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param q1 The first qubit.
	 * \param q2 The second qubit.
	 * \param angle Angle to rotate by.
	 */
	void npRotateX(unsigned q1, unsigned q2, Fp angle);

	/**
	 * \brief Perform an Y rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. This is equivalent to applying 
	 * \f$ e^{-i\theta (YX-XY)/2} \f$. 
	 *
	 * \ingroup qubits_gates

	 * \f[ 
	 *  npR_y(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -\sin(\theta/2) & 0 \\
	 *             0 & \sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param q1 The first qubit.
	 * \param q2 The second qubit.
	 * \param angle Angle to rotate by.
	 */
	void npRotateY(unsigned q1, unsigned q2, Fp angle);

	/**
	 * \brief Perform a Hadamard gate on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace.  
	 *
	 * \ingroup qubits_gates

	 * \f[ 
	 *  npH = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1/\sqrt{2} & 1/\sqrt{2} & 0 \\
	 *             0 & 1/\sqrt{2} & -1/\sqrt{2} & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * \param q1 The first qubit.
	 * \param q2 The second qubit.
	 */
	void npHadamard(unsigned q1, unsigned q2);

	
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
