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
 * \file time_multi.tpp
 * \brief Implementation of the multiple simulator version of Time class
 */

namespace qsl {

    template<typename Sim, Restrictions Res>
    template<typename R, typename S>
    void Time<Test::MultiSim, Sim, Res>::setMeta(
	const std::string & name, Results<R,S> & res)
    {
	res.addMeta("Testing '" + name + "' member function");
	res.addMeta("Using simulator '" + Sim::name + "'");
	res.addMeta("Number of qubits: " + std::to_string(nqubits));
	res.addMeta("Test length (function repeats per row): "
		    + std::to_string(test_len));
	res.addMeta("Using multiple random simulators "
		    "(one for each function evaluation)");
	res.addMeta("Format: space separated, one header row");
    };

    template<typename Sim, Restrictions Res>
    Time<Test::MultiSim, Sim, Res>::Time(unsigned nqubits_in,
					 std::size_t test_len_in,
					 std::size_t nsamples_in)
	: nqubits{nqubits_in}, test_len{test_len_in}, nsamples{nsamples_in}
    {
	// Print the test
	std::cout << "============================================"
		  << std::endl;

	std::cout << "Benchmarking " << Sim::name << std::endl;
	std::cout << "Generating " << test_len << " random state vectors"
		  << std::endl;

	// Make a list of random states and phases 
	for (std::size_t i = 0; i < test_len; i++) {
	    std::vector<complex<typename Sim::Fp_type>> temp
		= makeRandomState<typename Sim::Fp_type>(nqubits);
	    sim_list.push_back(temp);
	}

	std::cout << "Finished constructing benchmark class" << std::endl;
    }


    template<typename Sim, Restrictions Res>
    Results<unsigned, double>
    Time<Test::MultiSim, Sim, Res>::impl(const std::string & name,
					 OneQubitGate<Sim> fn)
    {		
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;
    
	Results<unsigned, double> res{{"qubit", name}};
	setMeta(name, res);
    
	for (unsigned n = 0; n < nqubits; n++) {
	    
	    qsl::Timer tmr;
	    
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(fn, sim_list[i], n);
	    }
	    tmr.stop();
	    double time = tmr.getElapsed();
	    
	    //printRow(std::to_string(n), time1, time2);
	    std::cout << n << ", " << time << std::endl;
	    res.addRow({n}, {time});
	}
	return res;
    }

    template<typename Sim, Restrictions Res>
    Results<unsigned, double>
    Time<Test::MultiSim, Sim, Res>::impl(const std::string & name,
					 OneQubitGate<Sim, typename Sim::Fp_type> fn)
    {	
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;

	Results<unsigned, double> res{{"qubit", name}};
	setMeta(name, res);
    	
	//std::vector<double> phase_list = makeRandomPhases(test_len);
	std::vector<typename Sim::Fp_type> phase_list
	    = makeRandomPhases<typename Sim::Fp_type>(test_len);
    
	for (unsigned n = 0; n < nqubits; n++) {
		      
	    qsl::Timer tmr;
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(fn, sim_list[i], n, phase_list[i]);
	    }
	    tmr.stop();
	    double time = tmr.getElapsed();
	    
	    //printRow(std::to_string(n), time1, time2);
	    std::cout << n << ", " << time << std::endl;
	    res.addRow({n}, {time});
	}
	return res;
    }

    template<typename Sim, Restrictions Res>
    Results<unsigned, double>
    Time<Test::MultiSim, Sim, Res>::impl(const std::string & name,
					 TwoQubitGate<Sim> fn)
    {	
	Results<unsigned, double> res{{"Control", "Target", name}};
	setMeta(name, res);

	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned m = n + 1; m < nqubits; m++) {
		qsl::Timer tmr;
		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(fn, sim_list[i], n, m);
		}
		tmr.stop();
		double time = tmr.getElapsed();
	
		//printRow(col1, time1, time2);
		std::cout << "(" << n << ", " << m << "), " << time << std::endl;
		res.addRow({n, m}, {time});
	    }
	}
	return res;
    }


    template<typename Sim, Restrictions Res>
    Results<unsigned, double>
    Time<Test::MultiSim, Sim, Res>::impl(
	const std::string & name,
	TwoQubitGate<Sim, typename Sim::Fp_type> fn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;

	Results<unsigned, double> res{{"Control", "Target", name}};
	setMeta(name, res);

	std::vector<typename Sim::Fp_type> phase_list
	    = makeRandomPhases<typename Sim::Fp_type>(test_len);

	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned m = n + 1; m < nqubits; m++) {

		qsl::Timer tmr;
	
		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(fn, sim_list[i], n, m, phase_list[i]);
		}
		tmr.stop();
		double time = tmr.getElapsed();
		
		std::cout << n << ", " << time << std::endl;
		res.addRow({n, m}, {time});	    
	    }
	}
	return res;
    }

    template<typename Sim, Restrictions Res>
    Results<std::size_t, double> Time<Test::MultiSim, Sim, Res>::sampleAll()
    {
	Results<std::size_t, double> r{{"nsamples", "sampleAll"}};
	setMeta("sampleAll", r);
    
	qsl::Timer tmr;

	// Work out the nsamples step size
	std::size_t step = nsamples/20;
    
	for (std::size_t n = 0; n < nsamples; n += step) {
    
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		sim_list[i].sampleAll(n);
	    }
	    tmr.stop();
	    double time = tmr.getElapsed();

	    r.addRow({n}, {time});
	    ///\todo Should be printing the row as well here (and in the others)
	}
    
	return r;
    }

}
