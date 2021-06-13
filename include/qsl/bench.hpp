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
 * \file bench.hpp
 * \brief Routines for testing the speed of the simulators
 * 
 */

#ifndef BENCH_HPP
#define BENCH_HPP

#include "utils/quantum.hpp"
#include "utils/timer.hpp"
#include <iostream>
#include <memory>
#include <functional>
#include <map>
#include "benchmark/results.hpp"

///\todo Write documentation
enum class Test
{
    SingleSim,
    MultiSim,
};

///\todo Get rid of all the Restrictions stuff. A much better
/// method for adding in number preserving states would be to
/// allow the state vector to be passed into the Time class
/// from the outside. This would make it the user's responsibility
/// to pick the right state vector.

template<Test T, typename Sim1, typename Sim2> class Compare;
template<Test T, typename Sim, Restrictions Res> class Time;

/**
 * \brief Run benchmark comparing Qubits gate and Quest gate
 *
 * \todo Move this compare to Compare<MultiSim>
 *
 * This benchmark compares the performance of a preselected gate
 * across all the qubits in two different simulator objects. First,  
 * a large number of random state vectors are generated, and then 
 * the gate is applied to all of the random state vectors in turn.
 * That operation is timed, and the results are printed. The function
 * repeats this for all the qubits (or all pairs of qubits for two
 * qubit gates).
 *
 * If the gate takes an argument (for example, rotateX), then a list
 * of random values is used for consecutive gate applications.
 *
 * The reason for using a new random state vector each time is to
 * prevent the compiler from making any simplifying assumptions
 * based on applying consecutive gates to the same state vector (for
 * example, two pauliX operations cancel out). It is not clear
 * whether the compiler could make use of any such properties. 
 *
 * The disadvantage of using multiple state vectors is that a new
 * state vector must be loaded from memory every time a get is
 * performed. This means that this benchmark may significantly
 * overestimate the runtime of a gate.
 *
 * \param fn A pointer to a gate member function of simulator F
 * \param gn A pointer to the same gate member function of simulator G
 * \param nqubits The number of qubits to include in state vectors
 * \param test_len The number of state vectors to generate and test
 * \param name A name to be printed along with the test results
 *
 * Note: the gate arguments are written like 
 * &Qubits<Type::Default>::pauliX
 *
 */
template<typename F, typename G, typename... ArgsF, typename... ArgsG>
void benchmark1(Gate<F, ArgsF...> fn, Gate<G, ArgsG...> gn,
		unsigned nqubits, int test_len);

/**
 * \brief Run benchmark comparing gates without random state vectors
 *
 * \todo Move this comment to Compare<SingleSim>
 *
 * This benchmark compares the performance of a preselected gate
 * across all the qubits in two different simulator objects. It differs
 * from bench1 in that only a single state vector is used for each
 * of the simulator objects, and the gates are applied to different
 * qubits in this single state vector.
 *
 * The use of a single state vector might reduce the chance of time-
 * consuming memory accesses which are constantly loading the different
 * state vectors from memory. In this benchmark, memory access times 
 * for each state vector should be less of a problem, so the runtime
 * should more accurately reflect the gate computation time.
 *
 * On the other hand, in simple cases (such as the pauliX gate), it
 * may be possible for the compiler to optimise away some of the 
 * operations if it is clever enough, which would invalidate the
 * benchmark.
 *
 * \param fn A pointer to a gate member function of simulator F
 * \param gn A pointer to the same gate member function of simulator G
 * \param nqubits The number of qubits to include in state vectors
 * \param test_len The number of repeated applications of the gate
 * \param name A name to be printed along with the test results
 *
 * Note: the gate arguments are written like 
 * &Qubits<Type::Default>::pauliX
 * 
 */
template<typename F, typename G, typename... ArgsF, typename... ArgsG>
void benchmark2(Gate<F, ArgsF...> fn, Gate<G, ArgsG...> gn,
		unsigned nqubits, int test_len);


/**
 * \brief Benchmark the measurement function
 *
 *
 */
template<typename F, typename G>
void benchmark1Measure(unsigned nqubits, int test_len);

/**
 * \brief Benchmark the measurement function
 *
 *
 */
template<typename F, typename G>
void benchmark2Measure(unsigned nqubits, int test_len);

/**
 * \brief Benchmkark the sampleAll function
 *
 * \todo Move this to Compare<MultiSim>::sampleAll
 *
 * The sampleAll function produces a list of measured outcomes from
 * the state of a simulation object, without collapsing the state 
 * vector. To benchmark the operation, this function generates
 * a large number of each kind of simulator object (test_len) and
 * then calls sampleAll(nsamples) on each object.
 *
 * Sometimes, if nsamples is large enough, test_len can be made
 * small because the sampling operating itself will take a long
 * time.
 *
 * \param nqubits The number of qubits to use for the tests
 * \param nsamples The number of samples to obtain from each sampleAll call
 * \param test_len The number of times to repeat each sampleAll call
 *
 */
template<typename F, typename G>
void benchmark1SampleAll(unsigned nqubits, int nsamples, int test_len = 1);

//#include "benchmark/bench1.tpp"
//#include "benchmark/bench2.tpp"

#endif

