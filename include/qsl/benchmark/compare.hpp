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
 * \file compare.hpp
 * \brief Put a brief description here...
 *
 * Put a detailed description here... 
 */

#ifndef QSL_COMPARE_HPP
#define QSL_COMPARE_HPP

#include "qsl/benchmark.hpp"

namespace qsl {
    /**
     * \brief Compare the speed of an operation between a
     * pair of randomly generated simulators.
     */
    template<typename Sim1, typename Sim2>
    class Compare<Test::SingleSim, Sim1, Sim2>
    {
	const unsigned nqubits;
	const std::size_t test_len; 
	const std::size_t nsamples; 

	Sim1 sim1;
	Sim2 sim2;

	/**
	 * \brief Set the metadata for the Results object
	 */
	template<typename R, typename S>
	void setMeta(const std::string & name, Results<R,S> & r);
    
	/**
	 * \brief Implementation for single qubit gates with no argument
	 */
	Results<unsigned, double>
	impl(const std::string & name, OneQubitGate<Sim1> fn,
	     OneQubitGate<Sim2> gn);
	/**
	 * \brief Implementation for single qubit gates with one double argument
	 */
	Results<unsigned, double>
	impl(const std::string & name,
	     OneQubitGate<Sim1, typename Sim1::Fp_type> fn,
	     OneQubitGate<Sim2, typename Sim2::Fp_type> gn);

	/**
	 * \brief Implementation for two qubit gates with no argument
	 */
	Results<unsigned, double>
	impl(const std::string & name, TwoQubitGate<Sim1> fn, TwoQubitGate<Sim2> gn);
    
	/**
	 * \brief Implementation for two qubit gates with one double argument
	 */
	Results<unsigned, double>
	impl(const std::string & name, TwoQubitGate<Sim1, typename Sim1::Fp_type> fn,
	     TwoQubitGate<Sim2, typename Sim2::Fp_type> gn);
    

    public:
    
	/**
	 * \brief Construct a Compare<Test::SingleSim> object
	 *
	 */
	Compare(unsigned nqubits, std::size_t test_len, std::size_t nsamples = 1);
    
	/**
	 * \brief Test pauliX
	 */
	Results<unsigned, double> pauliX() {
	    return impl("pauliX", &Sim1::pauliX, &Sim2::pauliX);
	}

	/**
	 * \brief Test rotateX
	 */
	Results<unsigned, double> rotateX() {
	    return impl("rotateX", &Sim1::rotateX, &Sim2::rotateX);
	}

	/**
	 * \brief Test phase
	 */
	Results<unsigned, double> phase() {
	    return impl("phase", &Sim1::phase, &Sim2::phase);
	}
    
	/**
	 * \brief Test controlNot
	 */
	Results<unsigned, double> controlNot() {
	    return impl("controlNot", &Sim1::controlNot,
			&Sim2::controlNot);
	}

	/**
	 * \brief Test controlPhase
	 */
	Results<unsigned, double> controlPhase() {
	    return impl("controlPhase", &Sim1::controlPhase,
			&Sim2::controlPhase);
	}

	Results<unsigned, double> measure();

	Results<std::size_t, double> sampleAll();
    
	///\todo Add sample, sampleAll, measureAll, prob
	///\todo Add the other gates
    };

    /**
     * \brief Compare the speed of an operation using multiple 
     * randomly generated simulators.
     */
    template<typename Sim1, typename Sim2>
    class Compare<Test::MultiSim, Sim1, Sim2>
    {
	const unsigned nqubits;
	const std::size_t test_len; 
	const std::size_t nsamples; 

	std::vector<Sim1> sim1_list;
	std::vector<Sim2> sim2_list;

	/**
	 * \brief Set the metadata for the Results object
	 */
	template<typename R, typename S>
	void setMeta(const std::string & name, Results<R,S> & r);
    
	/**
	 * \brief Implementation for single qubit gates with no argument
	 */
	Results<unsigned, double>
	impl(const std::string & name, OneQubitGate<Sim1> fn,
	     OneQubitGate<Sim2> gn);
	/**
	 * \brief Implementation for single qubit gates with one double argument
	 */
	Results<unsigned, double>
	impl(const std::string & name,
	     OneQubitGate<Sim1, typename Sim1::Fp_type> fn,
	     OneQubitGate<Sim2, typename Sim2::Fp_type> gn);

	/**
	 * \brief Implementation for two qubit gates with no argument
	 */
	Results<unsigned, double>
	impl(const std::string & name, TwoQubitGate<Sim1> fn, TwoQubitGate<Sim2> gn);
    
	/**
	 * \brief Implementation for two qubit gates with one double argument
	 */
	Results<unsigned, double>
	impl(const std::string & name, TwoQubitGate<Sim1, typename Sim1::Fp_type> fn,
	     TwoQubitGate<Sim2, typename Sim2::Fp_type> gn);
    
    public:
    
	/**
	 * \brief Construct a Compare<Test::MultiSim> object
	 *
	 */
	Compare(unsigned nqubits, std::size_t test_len, std::size_t nsamples = 1);
    
	/**
	 * \brief Test pauliX
	 */
	Results<unsigned, double> pauliX() {
	    return impl("pauliX", &Sim1::pauliX, &Sim2::pauliX);
	}

	/**
	 * \brief Test phase
	 */
	Results<unsigned, double> phase() {
	    return impl("phase", &Sim1::phase, &Sim2::phase);
	}
    
	/**
	 * \brief Test rotateX
	 */
	Results<unsigned, double> rotateX() {
	    return impl("rotateX", &Sim1::rotateX, &Sim2::rotateX);
	}

	/**
	 * \brief Test controlNot
	 */
	Results<unsigned, double> controlNot() {
	    return impl("controlNot", &Sim1::controlNot,
			&Sim2::controlNot);
	}

	/**
	 * \brief Test controlPhase
	 */
	Results<unsigned, double> controlPhase() {
	    return impl("controlPhase", &Sim1::controlPhase,
			&Sim2::controlPhase);
	}

	Results<unsigned, double> measure();

	Results<std::size_t, double> sampleAll();
    
	///\todo Add sample, sampleAll, measureAll, prob
	///\todo Add the other gates
    };
}
#include "qsl/benchmark/compare_multi.tpp"
#include "qsl/benchmark/compare_single.tpp"

#endif 
