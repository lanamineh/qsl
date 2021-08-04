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
 * \file time_single.tpp
 * \brief Implementation of the single simulator version of Time class
 */

namespace qsl {

    template<typename Sim, Restrictions Res>
    template<typename R, typename S>
    void Time<Test::SingleSim, Sim, Res>::setMeta(
	const std::string & name, Results<R,S> & res)
    {
	res.addMeta("Testing '" + name + "' member function");
	res.addMeta("Using simulator '" + Sim::name + "'");
	res.addMeta("Number of qubits: " + std::to_string(nqubits));
	res.addMeta("Test length (function repeats per row): "
		    + std::to_string(test_len));
	res.addMeta("Using a single random simulator for every function evaluation");
	res.addMeta("Format: space separated, one header row");
    };

    template<typename Sim, Restrictions Res>
    Time<Test::SingleSim, Sim, Res>::Time(unsigned nqubits_in,
					  std::size_t test_len_in,
					  std::size_t nsamples_in)
	: nqubits{nqubits_in}, test_len{test_len_in}, nsamples{nsamples_in},
	  sim{nqubits}
    {
	// Print the test
	std::cout << "============================================"
		  << std::endl;

	std::cout << "Benchmarking " << Sim::name << std::endl;
    
	// Make a random state vector to start with
    
	std::vector<complex<typename Sim::Fp_type>> temp
	    = RandomStateGen<typename Sim::Fp_type, Res>(nqubits).get();
	sim.setState(temp);

	std::cout << "Finished constructing benchmark class" << std::endl;
    }

    template<typename Sim, Restrictions Res>
    Results<unsigned, double> Time<Test::SingleSim, Sim, Res>::impl(
	const std::string & name,
	OneQubitGate<Sim> fn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;
	
	using R = Results<unsigned, double>;
	R r{{"qubit", name}};
	setMeta(name, r);
    
	for (unsigned n = 0; n < nqubits; n++) {
		      
	    qsl::Timer tmr;
		      
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(fn, sim, n);
	    }
	    tmr.stop();
	    double time = tmr.getElapsed();

	    //printRow(std::to_string(n), time1, time2);
	    std::cout << n << ", " << time << std::endl;
	    r.addRow({n}, {time});

	}

	return r;
    }

    template<typename Sim, Restrictions Res>
    Results<unsigned, double> Time<Test::SingleSim, Sim, Res>::impl(
	const std::string & name,
	OneQubitGate<Sim, typename Sim::Fp_type> fn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;
    
	using R = Results<unsigned, double>;
	R r{{"qubit", name}};
	setMeta(name, r);
    
	std::vector<typename Sim::Fp_type> phase_list
	    = makeRandomPhases<typename Sim::Fp_type>(nqubits);
	
	for (unsigned n = 0; n < nqubits; n++) {
		      
	    qsl::Timer tmr;
		      
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(fn, sim, n, phase_list[n]);
	    }
	    tmr.stop();
	    double time = tmr.getElapsed();
	    
	    //printRow(std::to_string(n), time1, time2);
	    std::cout << n << ", " << time << std::endl;	
	    r.addRow({n}, {time});
	}

	return r;
    }

    template<typename Sim, Restrictions Res>
    Results<unsigned, double> Time<Test::SingleSim, Sim, Res>::impl(
	const std::string & name,
	TwoQubitGate<Sim> fn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;
    
	using R = Results<unsigned, double>;
	R r{{"Control", "Target", name}};
	setMeta(name, r);
	
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned m = n + 1; m < nqubits; m++) {

		qsl::Timer tmr;
	
		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(fn, sim, n, m);
		}
		tmr.stop();
		double time = tmr.getElapsed();
		
		std::string col1 = "(" + std::to_string(n) + ", ";
		col1 += std::to_string(m) + ")";
		//printRow(col1, time1, time2);
		std::cout << "(" << n << ", " << m << "), " << time << std::endl;
		r.addRow({n, m}, {time});
	    }
	}

	return r;
    }

    template<typename Sim, Restrictions Res>
    Results<unsigned, double> Time<Test::SingleSim, Sim, Res>::impl(
	const std::string & name,
	TwoQubitGate<Sim, typename Sim::Fp_type> fn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;
    
	using R = Results<unsigned, double>;
	R r{{"Control", "Target", name}};
	setMeta(name, r);

	std::vector<typename Sim::Fp_type> phase_list
	    = makeRandomPhases<typename Sim::Fp_type>(nqubits*nqubits);
    
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned m = n + 1; m < nqubits; m++) {
	    
		qsl::Timer tmr;
	
		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(fn, sim, n, m, phase_list[nqubits*n + m]);
		}
		tmr.stop();
		double time = tmr.getElapsed();

		std::string col1 = "(" + std::to_string(n) + ", ";
		col1 += std::to_string(m) + ")";
		//printRow(col1, time1, time2);
		std::cout << "(" << n << ", " << m << "), " << time << std::endl;
		r.addRow({n, m}, {time});
	    }
	}
	return r;
    }

    template<typename Sim, Restrictions Res>
    Results<std::size_t, double> Time<Test::SingleSim, Sim, Res>::sampleAll()
    {
	Results<std::size_t, double> r{{"nsamples", "sampleAll"}};
	setMeta("sampleAll", r);
    
	qsl::Timer tmr;

	// Work out the nsamples step size
	std::size_t step = nsamples/20;
    
	for (std::size_t n = 0; n < nsamples; n += step) {
    
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		sim.sampleAll(n);
	    }
	    tmr.stop();
	    double time = tmr.getElapsed();

	    r.addRow({n}, {time});
	    ///\todo Should be printing the row as well here (and in the others)
	}
    
	return r;
    }

}
