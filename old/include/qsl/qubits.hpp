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
 * \file qubits.hpp
 * \brief Contains the main simulator object for manipulating qubits
 * 
 */

#ifndef QSL_QUBITS_HPP
#define QSL_QUBITS_HPP

#include <concepts>

///\todo At the moment this is only required because of the use
/// of fubiniStudy in the template function below. Maybe the
/// fubiniStudy function should be moved somewhere else.
#include <qsl/utils.hpp>

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
	OmpNP, ///< Number preserving simulator with multithreading
	Resize, ///< Resizeable simulator
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
    template<Type T = Type::Default, std::floating_point Fp = double> class Qubits;

}
    
///\todo Find a good structure for including these files
#include "qubits/default.hpp"
#include "qubits/omp.hpp"
#include "qubits/np.hpp"
#include "qubits/omp_np.hpp"
#include "qubits/resize.hpp"

#include "concepts.hpp"

namespace qsl {

    ///\todo Maybe these typedefs are a bad idea?
    template<typename Sim, typename... Args>
    using Gate = void(Sim::*)(Args...);

    template<typename Sim, typename... Args>
    using OneQubitGate = Gate<Sim, unsigned, Args...>;

    template<typename Sim, typename... Args>
    using TwoQubitGate = Gate<Sim, unsigned, unsigned, Args...>;

    ///\todo Maybe this should go in utils?
    template<Simulator S1, Simulator S2>
    requires std::is_same_v<typename S1::Fp_type, typename S2::Fp_type>
    S1::Fp_type fubiniStudy(const S1 & s1, const S2 & s2)
    {
	return qsl::fubiniStudy(s1.getState(), s2.getState());
    }
    
}
#endif
