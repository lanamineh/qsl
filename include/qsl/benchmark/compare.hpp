/* 
 * Copyright (C) 2020 Lana Mineh and John Scott.
 *
 * This file is part of %PROJNAME%, the quantum computer simulator.
 *
 * %PROJNAME% is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * %PROJNAME% is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with %PROJNAME%.  If not, see <https://www.gnu.org/licenses/>.
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
