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
 * \file verify.hpp
 * \brief Testing whether the simulators agree
 * 
 * It is important to check that all the simulator objects give
 * the same results when gates are applied and measurements are
 * performed. For gates and other deterministic operations, the 
 * check is simple: whether the state vectors agree before and 
 * after the gate is applied. For probabilistic operations, such
 * as measurement, it is necessary to test whether the simulators
 * give the same average results after a large number of 
 * measurements have been performed.
 *
 * This file contains a class, Verify, for performing these 
 * comparisons. Any pair of simulator objects may be tested
 * against each other, so that it is possible to cross 
 * reference all the simulators and check that they agree.
 */

#ifndef VERIFY_HPP
#define VERIFY_HPP

#include "qsl/concepts.hpp"
#include "qsl/utils/quantum.hpp"
#include "qsl/utils/misc.hpp"
#include "qsl/verify/checkers.hpp"
#include "qsl/verify/stategen.hpp"

namespace qsl {

    /**
     * \brief Verify two simulators against one another
     *
     * This class can be used to check that two simulators give the same
     * results across a range of different functions. It takes several
     * template parameters, which control what it does. The first two
     * parameters specify the simulators under test. The next parameter
     * is a class which controls how the simulators are initialised. The
     * other template parameters are classes which test different aspects
     * of the simulators.
     *
     * Calling the checkAll method will pass through the list of checkers,
     * re-initialising the state using the StateGen class, and checking the
     * simulators using the current Checker class in the Checkers... parameter.
     * The function returns after everything has been checked.
     *
     * It is possible to customise the behaviour of the Verify class by writing
     * new state generators and checkers. 
     *
     */
    template<qsl::Simulator Sim1, qsl::Simulator Sim2, qsl::StateGenerator Gen,
	     template<typename,typename> class... Checkers>
    // Constrain the Checkers type using the SimChecker concept
    ///\todo This causes a bug at the moment in the checkers
    //requires (qsl::SimChecker<Checkers<Sim1,Sim2>,Sim1,Sim2>>...)
    class Verify
    {
	Gen init;
	std::tuple<Checkers<Sim1, Sim2>...> checkers;

	/**
	 * \brief Initialise two simulators and run the checker
	 */
	template<template<typename,typename> class Checker>
	Checker<Sim1,Sim2>::ResultData initAndCheck(std::ostream & os)
	    {
		os << "Resetting state vector" << std::endl;
		std::unique_ptr<Sim1> sim1
		    = std::make_unique<Sim1>(init.getState());
		std::unique_ptr<Sim2> sim2
		    = std::make_unique<Sim2>(init.getState());
		std::get<Checker<Sim1,Sim2>>(checkers).bind(sim1, sim2);
		return std::get<Checker<Sim1,Sim2>>(checkers).checkAll(os);
	    }
    
    public:

	template<typename... Args>
	void configureState(Args... args)
	    {
		init.configureState(args...);
	    }

	template<template<typename,typename> class Checker, typename... Args>
	void configureChecker(Args... args)
	    {
		std::get<Checker<Sim1,Sim2>>(checkers).configureChecker(args...);
	    }
	
	/**
	 * \brief Run one of the checkers and get results
	 *
	 * This function is for the unit test, where you want to automatically
	 * check the results of the checker rather than print them out.
	 *
	 * The returned type is specific to the checker and contains the
	 * results of the tests. Look at the ResultData struct inside the
	 * checker.
	 */
	template<template<typename,typename> class Checker, typename... Args>
	Checker<Sim1,Sim2>::ResultData check()
	    {
		qsl::NullStream nul;
		return initAndCheck<Checker>(nul); 
	    }

	/** 
	 * \brief Use this function to print the results to the screen
	 *
	 */
	void checkAll()
	    {
		// Execute all the checks (C++17 fold expression)
		(initAndCheck<Checkers>(std::cout), ...); 
	    }	  
    };

}
#endif
