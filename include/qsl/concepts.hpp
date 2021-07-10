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
 * \file concepts.hpp
 * \brief Concepts for different simulator types
 * 
 */

#ifndef CONCEPTS_HPP
#define CONCEPTS_HPP

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
	sim.phase(targ,param);
	sim.controlPhase(ctrl,targ,param);
	sim.swap(ctrl,targ); // Actually q1 and q2, but still unsigned
	sim.controlZ(ctrl,targ);
	sim.rotateZ(targ,param);
    };
    
    template<typename Sim>
    concept HasNonNPGates = requires (Sim sim, unsigned ctrl, unsigned targ,
				      Sim::Fp_type param) {
	sim.hadamard(targ);
	sim.pauliX(targ);
	sim.rotateX(targ,param);
	sim.controlNot(ctrl,targ);
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
