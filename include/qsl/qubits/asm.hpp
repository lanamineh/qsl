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
 * \file qubits/asm.hpp
 * \brief Contains the assembly written simulator object for manipulating qubits
 * 
 */

#ifndef QUBITS_ASM_HPP
#define QUBITS_ASM_HPP

/**
 * \defgroup gates Quantum gates
 * \brief One- and two- qubit gates in the Qubits class
 */

/**
 * \defgroup meas Measurement
 * \brief Functions related to measurement, probabilities, etc.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"

/**
 * \brief Assembly quantum simulator object
 *
 * This version of the Qubits object uses gates and other
 * methods that are written using x86 assembly. It is
 * predominantly useful as an optimising tool for the 
 * other classes. Threading is not used.
 *
 */
template<>
class Qubits<Type::Asm, double>
{
    const unsigned nqubits; 
    const std::size_t dim; ///< The length of the state vector
    std::vector<complex<double>> state; ///< State vector for the qubits
    
    qsl::Random<double> random;

    struct Dist {
	std::size_t index;
	double prob;
    };
    
    /// Helper function for measurement related functions
    std::size_t drawSample(const std::vector<Dist> & dist);
    
public:

    const static std::string name; ///< The name of this simulator
    
    /// Expose Fp to users of the class
    using Fp_type = double; /// \todo Is there a less bad way to do this?
    
    /// Initialise the class with a fixed number of qubits
    Qubits(unsigned nqubits);

    /// Initialise the class from a pre-prepared state vector
    Qubits(const std::vector<complex<Fp_type>> & state);

    /// Copy constructor
    Qubits(const Qubits & ) = default;
    
    /// Copy assignment operator
    void operator = (const Qubits & old);

    /// Return the number of qubits
    unsigned getNumQubits() const;
    
    /// Reset to the all-zero state
    void reset();

    /// Print the state vector
    void print() const;

    /// Get the state vector associated to the qubits
    std::vector<complex<Fp_type>> getState() const;

    /// Set the state vector (i.e. re-initialise the state vector)
    void setState(const std::vector<complex<Fp_type>> & state);

    /// Set the state vector to a computational basis state
    void setBasisState(std::size_t index);
    
    /**
     * \brief Rotate around the x-axis of the Bloch sphere \f$ e^{-i\theta X/2} \f$
     *
     * \ingroup gates
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
    void rotateX(unsigned targ, Fp_type angle);

    /**
     * \brief Apply the Pauli X gate to qubit number targ.
     *
     * \ingroup gates
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
     * \brief Measure a qubit and collapse the state to its outcome.
     *
     * \ingroup meas
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
     * \ingroup meas
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
     * \ingroup meas
     *
     * \param targ The qubit to calculate the probability for.
     * \param outcome The outcome (0 or 1) we are calculating the probability of.
     * \return The probability of the qubit being measured in the given outcome.
     */
    double prob(unsigned targ, unsigned outcome) const;

    /**
     * \brief Perform a post-selection measurement. The state is collapsed 
     * to the given outcome for the given qubit. 
     *
     * \ingroup meas
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
     * \ingroup meas
     *
     * \param targ The qubit to measure.
     * \param nsamples The number of samples to draw.
     * \return A vector containing all the measurement results.
     */
    std::vector<int> sample(unsigned targ, std::size_t nsamples);
    
    /**
     * \brief Sample measurement outcome for all of the qubits 
     * at once multiple times.
     *
     * \ingroup meas
     *
     * Measuring all of the qubits at once is the same as measuring them 
     * one by one. This function implements a very efficient way of
     * measuring all the qubits multiple times and should be used 
     * instead of the measure function where possible. Note that this 
     * function does not modify the state vector.
     * 
     * \param nmeas The number of measurements to perform.
     * \return A vector containing all of the measurement results, each
     *         result is returned as the computational basis state index. 
     */
    std::vector<std::size_t> sampleAll(std::size_t nsamples);
};

#endif
