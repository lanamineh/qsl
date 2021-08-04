/*
 *  Copyright 2021 Lana Mineh and John Scott
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
 * \file concepts.hpp
 * \brief Concepts for different simulator types
 * 
 */

#ifndef QSL_CONCEPTS_HPP
#define QSL_CONCEPTS_HPP

#include <concepts>
#include <iostream>

namespace qsl {

    /**
     * \brief Simulator object concept
     *
     * All simulators should conform to this specification.
     * A simulator object T must:
     *
     * 1) have a constructor T(std::size_t) that constructs a 
     *    default object from a given number of qubits;
     * 2) have a constructor T(std::vector<qsl::complex<Fp>>)
     *    which constructs a simulator object from a given
     *    state vector;
     * 3) be copy-constructible and copy-assignable;
     * 4) have print, reset, getNumQubits and getState methods;
     * 5) have a setState method.
     *
     * 
     */
    template<typename T>
    concept Simulator = std::is_floating_point<typename T::Fp_type>::value
	&& std::is_constructible<std::size_t>::value
	&& std::is_constructible<std::vector<
				     qsl::complex<typename T::Fp_type>>>::value
	&& std::is_copy_constructible<T>::value
	&& std::is_copy_assignable<T>::value
    	&& requires(T t, std::vector<qsl::complex<typename T::Fp_type>> state)
    {
    	// Required member functions
    	t.reset();
    	t.getNumQubits();
    	t.print();
    	t.getState();
	t.setState(state);
    };
    
    template<typename Sim>
    concept HasNPGates = requires (Sim sim, unsigned ctrl, unsigned targ,
				   Sim::Fp_type param) {
	// One qubit gates
	sim.pauliZ(targ); 
	sim.rotateZ(targ,param);
	sim.phase(targ,param);	

	// Controlled two qubit gates
	sim.controlZ(ctrl,targ);
	sim.controlRotateZ(ctrl,targ,param);
	sim.controlPhase(ctrl,targ,param);

	// Number preserved two qubit gates
	sim.swap(ctrl,targ); // Actually q1 and q2, but still unsigned
	sim.fswap(ctrl,targ);
	sim.npRotateX(ctrl,targ,param);
	sim.npRotateY(ctrl,targ,param);
	sim.npHadamard(ctrl,targ);
    };
    
    template<typename Sim>
    concept HasNonNPGates = requires (Sim sim, unsigned ctrl, unsigned targ,
				      Sim::Fp_type param) {

	// One qubit gates
	sim.pauliX(targ);
	sim.pauliY(targ);
	sim.rotateX(targ,param);
	sim.rotateY(targ,param);
	sim.hadamard(targ);

	// Controlled two qubit gates
	sim.controlNot(ctrl,targ);
	sim.controlY(ctrl,targ);
	sim.controlRotateX(ctrl,targ,param);	
	sim.controlRotateY(ctrl,targ,param);
	sim.controlHadamard(ctrl,targ);
    };

    template<typename Sim>
    concept HasMeasurement = requires (Sim sim, unsigned targ,
				       unsigned outcome, std::size_t nsamples) {
	///\todo Change some of the returned types to use Fp_type
	{ sim.measure(targ) } -> std::same_as<int>;
	{ sim.measureAll() } -> std::same_as<std::size_t>;
	{ sim.prob(targ,outcome) } -> std::convertible_to<double>;
	{ sim.postselect(targ,outcome) } -> std::convertible_to<double>;
	{ sim.sample(targ,nsamples) };// -> std::same_as<std::vector<std::size_t>>;
	{ sim.sampleAll(nsamples) };// -> std::same_as<std::map<std::size_t,std::size_t>>;
    };

    /// Definition of a type 
    template<typename Sim>
    concept HasAllGates = HasNPGates<Sim> && HasNonNPGates<Sim>;

    template<typename Sim>
    concept NPSimulator = HasMeasurement<Sim> && HasNPGates<Sim>;

    /// Just in case you want...
    template<typename Sim>
    concept NonNPSimulator = HasMeasurement<Sim> && HasNonNPGates<Sim>;

}
#endif
