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
 * \file thread.hpp
 * \brief Contains the threaded implementation of the Qubits simulator
 *
 */

#ifndef THREAD_HPP
#define THREAD_HPP

/**
 * \defgroup qubits_thread_constructors Constructors
 * \brief Constructors for the threaded simulator object
 */

/**
 * \defgroup qubits_thread_gates Gates
 * \brief Gates for the threaded simulator object
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include "qsl/utils/complex.hpp"
#include "qsl/utils/random.hpp"
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <iostream>

#include "qsl/threads/pool.hpp"

template<typename Fp>
class Qubits<Type::Thread, Fp> {

    const unsigned nqubits; 
    const std::size_t dim; ///< The length of the state vector
    std::vector<complex<Fp>> state; ///< State vector for the qubits
    
    Pool pool; ///< Thread pool

    qsl::Random<Fp> random;
        
    struct Dist {
	std::size_t index;
	double prob;
    };

    void safePrint(const std::string & str) {
	static std::mutex cout_mutex;       
	std::lock_guard<std::mutex> lock(cout_mutex);
	std::cout << str << std::endl;
    }

    void pauliX_inner(std::size_t start, std::size_t end, std::size_t k);
    void pauliX_outer(std::size_t start, std::size_t end, std::size_t k);
    void pauliX_0(std::size_t start, std::size_t end);
    void pauliX_1(std::size_t start, std::size_t end);
    void pauliX_2(std::size_t start, std::size_t end);
    void pauliX_3(std::size_t start, std::size_t end);
	
public:

    const static std::string name; ///< The name of this simulator
    
    /// Expose Fp to users of the class
    using Fp_type = Fp; ///< \todo Is there a less bad way to do this?
    
    /**
     * \brief Initialise the class with a fixed number of qubits
     * \ingroup qubits_thread_constructors
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
     * \ingroup qubits_thread_constructors
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
     * \ingroup qubits_thread_constructors
     *
     * You can make copies of this object by constructing from an 
     * object of the same type.
     */
    Qubits(const Qubits & );

    /**
     * \brief Copy-assignment operator
     * \ingroup qubits_thread_constructors
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
     * \brief Apply the Pauli X gate to qubit number targ.
     *
     * \ingroup qubits_thread_gates
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
     * \brief Perform the controlled Not (CNOT) gate on two qubits.
     *
     * \ingroup qubits_thread_gates
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


};

#endif 
