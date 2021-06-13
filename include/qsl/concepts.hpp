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

template<typename Sim>
concept HasNPGates = requires (Sim sim, unsigned ctrl, unsigned targ,
			       Sim::Fp_type param) {
    sim.phase(targ,param);
    sim.controlPhase(ctrl,targ,param);
    // Swap too
};

template<typename Sim>
concept HasNonNPGates = requires (Sim sim, unsigned ctrl, unsigned targ,
				  Sim::Fp_type param) {
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

#endif
