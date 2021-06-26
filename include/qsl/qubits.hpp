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
 * \file qubits.hpp
 * \brief Contains the main simulator object for manipulating qubits
 * 
 */

#ifndef QUBITS_HPP
#define QUBITS_HPP

namespace qsl {
/**
 * \defgroup gates Quantum gates
 * \brief One- and two- qubit gates in the Qubits class
 */

/**
 * \defgroup meas Measurement
 * \brief Functions related to measurement, probabilities, etc.
 */

/**
 * \brief An enum for distinguishing between different simulator 
 * implementations
 *
 * A value from this enum should be passed as the first template parameter
 * to the Qubits class to determine which simulator should be used.
 */
enum class Type {
    Default, ///< The default simulator object (no inline assembly or threading)
    Omp, ///< Simulator supporting multithreading based on OpenMP
    NP, ///< Number preserving simulator
};

/**
 * \brief The main quantum simulator object
 *
 * There are several different types of simulator
 * objects, depending on the implementation. The
 * type of simulator is chosen by specifying the
 * implementation using the Type enum, and specifying
 * the simulator precision Fp. 
 * 
 */
template<Type T = Type::Default, typename Fp = double> class Qubits;

}
    
///\todo Find a good structure for including these files
#include "qubits/default.hpp"
#include "qubits/omp.hpp"
#include "qubits/np.hpp"

namespace qsl {

///\todo Maybe these typedefs are a bad idea?
template<typename Sim, typename... Args>
using Gate = void(Sim::*)(Args...);

template<typename Sim, typename... Args>
using OneQubitGate = Gate<Sim, unsigned, Args...>;

template<typename Sim, typename... Args>
using TwoQubitGate = Gate<Sim, unsigned, unsigned, Args...>;

}
#endif
