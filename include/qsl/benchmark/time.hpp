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
 * \file time.hpp
 * \brief Put a brief description here...
 *
 * Put a detailed description here... 
 */

#ifndef TIME_HPP
#define TIME_HPP

#include "qsl/bench.hpp"

enum class Function
{
    pauliX,
    rotateX,
    phase,
    controlPhase,
    controlNot,
};

/**
 * \brief Run the timing test
 */
template<Function Fn>
struct FnData
{
    static const unsigned col_num; ///< The column number with timing results
    static const std::string name; ///< The function name
    FnData(); ///< Constructor

    template<typename Sim, Restrictions Res>
    static Results<unsigned,double> run(Time<Test::SingleSim,Sim,Res> & time);
};

// Pauli-X
template<> const unsigned FnData<Function::pauliX>::col_num = 1;
template<> const std::string FnData<Function::pauliX>::name = "pauliX";
template<> template<typename Sim, Restrictions Res>
Results<unsigned,double>
FnData<Function::pauliX>::run(Time<Test::SingleSim,Sim,Res> & time)
{
    return time.pauliX();
}

// Phase
template<> const unsigned FnData<Function::phase>::col_num = 1;
template<> const std::string FnData<Function::phase>::name = "phase";
template<> template<typename Sim, Restrictions Res>
Results<unsigned,double>
FnData<Function::phase>::run(Time<Test::SingleSim,Sim,Res> & time)
{
    return time.phase();
}


// Rotate-X
template<> const unsigned FnData<Function::rotateX>::col_num = 1;
template<> const std::string FnData<Function::rotateX>::name = "rotateX";
template<> template<typename Sim, Restrictions Res>
Results<unsigned,double>
FnData<Function::rotateX>::run(Time<Test::SingleSim,Sim,Res> & time)
{
    return time.rotateX();
}


// Control NOT
template<> const unsigned FnData<Function::controlNot>::col_num = 2;
template<> const std::string FnData<Function::controlNot>::name = "controlNot";
template<> template<typename Sim, Restrictions Res>
Results<unsigned,double>
FnData<Function::controlNot>::run(Time<Test::SingleSim,Sim,Res> & time)
{
    return time.controlNot();
}

// Control Phase
template<> const unsigned FnData<Function::controlPhase>::col_num = 2;
template<> const std::string FnData<Function::controlPhase>::name = "controlPhase";
template<> template<typename Sim, Restrictions Res>
Results<unsigned,double>
FnData<Function::controlPhase>::run(Time<Test::SingleSim,Sim,Res> & time)
{
    return time.controlPhase();
}

/**
 * \brief Timed sweep over number of qubits
 *
 * This function times each simulator in the template argument list, for
 * each number of qubits up to a limit specified in the argument list. The
 * function under test is specified in the argument list, as is a name for
 * the output file. 
 *
 * The function creates a Time<SingleSim> object for each total number of 
 * qubits n up to nqubits_max. From this object, a Results class is obtained
 * containing the time that the functino fn takes on each qubit in the state
 * vector. These results are averaged and stored in a single row in the output
 * file. This row represents the average time it takes to perform the function
 * on any qubit (or pair of qubits) for a given total number of qubits n.
 *
 *
 */
template<Function Fn, Restrictions Res, typename... Sims>
void nqubitSweep(std::string results_file, unsigned len_factor = 20)
{
    unsigned nqubits_max = 19;
    Results<unsigned, double> results{ {"nqubits", Sims::name... } };
    const std::size_t base_len = (std::size_t)1 << len_factor;

    // Loop over total number of qubits
    for (unsigned n = 2; n < nqubits_max; n++) {

	std::size_t test_len = std::max((std::size_t)1,
					base_len / ((std::size_t)1 << n)); 
	std::cout << "Test len = " << test_len << std::endl;

	// Make a tuple of all the timing objects for each simulator
	std::tuple<Time<Test::SingleSim, Sims, Res>...> time{
	    Time<Test::SingleSim, Sims, Res>{n,test_len}...
	};
	
	// Need to choose the right column for one or two qubit gates
	std::vector<Results<unsigned,double>> timing_results{
	    FnData<Fn>::run(std::get<Time<Test::SingleSim, Sims, Res>>(time))...
	};
	
	// Average the time over each row in the results
	std::vector<double> averages;
	for (std::size_t i = 0; i < timing_results.size(); i++) {
	    unsigned col_num = FnData<Fn>::col_num;
	    averages.push_back(timing_results[i].mean<double>(col_num)/test_len);
	}	
	results.addRow({n}, averages);
    }

    results.addMeta("Averaged gate times for " + FnData<Fn>::name);
    results.addMeta("Results are in seconds");
    
    results.writeToFile(results_file);
    std::cout << std::endl << "Written results to " << results_file << std::endl;
}


/**
 * \brief Time an operation using a single random simulator,
 */
template<typename Sim, Restrictions Res>
class Time<Test::SingleSim, Sim, Res>
{
    const unsigned nqubits;
    const std::size_t test_len; 
    const std::size_t nsamples; 

    Sim sim;
    
    /**
     * \brief Set the metadata for the Results object
     */
    template<typename R, typename S>
    void setMeta(const std::string & name, Results<R,S> & r);
    
    /**
     * \brief Implementation for single qubit gates with no argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   OneQubitGate<Sim> fn);

    /**
     * \brief Implementation for single qubit gates with one double argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   OneQubitGate<Sim, typename Sim::Fp_type> fn);

    /**
     * \brief Implementation for two qubit gates with no argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   TwoQubitGate<Sim> fn);
    
    /**
     * \brief Implementation for two qubit gates with one double argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   TwoQubitGate<Sim, typename Sim::Fp_type> fn);
    

public:
    
    /**
     * \brief Construct a Time object
     *
     */
    Time(unsigned nqubits, std::size_t test_len, std::size_t nsamples = 1);
    
    /// Print
    int print() const {
	std::cout << "Time (SingleSim)" << std::endl;
	std::cout << "- Number of qubits = " << nqubits << std::endl;
	std::cout << "- Test length = " << test_len << std::endl;
	std::cout << "- Number of samples = " << nsamples << std::endl;
	return 0;
    }

    /**
     * \brief Test pauliX
     */
    Results<unsigned, double> pauliX() {
	return impl("pauliX", &Sim::pauliX);
    }

    /**
     * \brief Test rotateX
     */
    Results<unsigned, double> rotateX() {
	return impl("rotateX", &Sim::rotateX);
    }

    /**
     * \brief Test phase
     */
    Results<unsigned, double> phase() {
	return impl("phase", &Sim::phase);
    }
    
    /**
     * \brief Test controlNot
     */
    Results<unsigned, double> controlNot() {
	return impl("controlNot", &Sim::controlNot);
    }

    /**
     * \brief Test controlPhase
     */
    Results<unsigned, double> controlPhase() {
	return impl("controlPhase", &Sim::controlPhase);
    }

    Results<std::size_t, double> sampleAll();

};

/**
 * \brief Time an operation using multiple randomly generated simulators
 */
template<typename Sim, Restrictions Res>
class Time<Test::MultiSim, Sim, Res>
{
    const unsigned nqubits;
    const std::size_t test_len; 
    const std::size_t nsamples; 

    std::vector<Sim> sim_list;

    /**
     * \brief Set the metadata for the Results object
     */
    template<typename R, typename S>
    void setMeta(const std::string & name, Results<R,S> & r);
    
    /**
     * \brief Implementation for single qubit gates with no argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   OneQubitGate<Sim> fn);

    /**
     * \brief Implementation for single qubit gates with one double argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   OneQubitGate<Sim, typename Sim::Fp_type> fn);

    /**
     * \brief Implementation for two qubit gates with no argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   TwoQubitGate<Sim> fn);
    
    /**
     * \brief Implementation for two qubit gates with one double argument
     */
    Results<unsigned, double> impl(const std::string & name,
				   TwoQubitGate<Sim, typename Sim::Fp_type> fn);

    
public:
    
    /**
     * \brief Construct a Compare<Test::MultiSim> object
     *
     */
    Time(unsigned nqubits, std::size_t test_len, std::size_t nsamples = 1);
    
    /**
     * \brief Test pauliX
     */
    Results<unsigned, double> pauliX() {
	return impl("pauliX", &Sim::pauliX);
    }

    /**
     * \brief Test phase
     */
    Results<unsigned, double> phase() {
	return impl("phase", &Sim::phase);
    }
    
    /**
     * \brief Test rotateX
     */
    Results<unsigned, double> rotateX() {
	return impl("rotateX", &Sim::rotateX);
    }

    /**
     * \brief Test controlNot
     */
    Results<unsigned, double> controlNot() {
	return impl("controlNot", &Sim::controlNot);
    }

    /**
     * \brief Test controlPhase
     */
    Results<unsigned, double> controlPhase() {
	return impl("controlPhase", &Sim::controlPhase);
    }

    Results<std::size_t, double> sampleAll();

};

#include "qsl/benchmark/time_multi.tpp"
#include "qsl/benchmark/time_single.tpp"

#endif 
