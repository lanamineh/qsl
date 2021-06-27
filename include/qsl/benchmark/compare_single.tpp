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
 * \file compare_single.tpp
 * \brief Implementation of the single simulator version of Compare class
 */

namespace qsl {

    template<typename Sim1, typename Sim2>
    template<typename R, typename S>
    void Compare<Test::SingleSim, Sim1, Sim2>::setMeta(
	const std::string & name, Results<R,S> & res)
    {
	res.addMeta("Testing '" + name + "' member function");
	res.addMeta("Comparing '" + Sim1::name + "' and '" + Sim2::name + "'");
	res.addMeta("Number of qubits: " + std::to_string(nqubits));
	res.addMeta("Test length (function repeats per row): "
		    + std::to_string(test_len));
	res.addMeta("Using a single random simulator for every function evaluation");
	res.addMeta("Format: space separated, one header row");
    };

    template<typename Sim1, typename Sim2>
    Compare<Test::SingleSim, Sim1, Sim2>::Compare(unsigned nqubits_in,
						  std::size_t test_len_in,
						  std::size_t nsamples_in)
	: nqubits{nqubits_in}, test_len{test_len_in}, nsamples{nsamples_in},
	  sim1{nqubits}, sim2{nqubits}
    {
	// Print the test
	std::cout << "============================================"
		  << std::endl;

	std::cout << "Benchmarking " << std::endl
		  << Sim1::name << " (1) and " << Sim2::name
		  << " (2)." << std::endl;
    
	// Make a random state vector to start with
	std::vector<complex<typename Sim1::Fp_type>> temp
	    = makeRandomState<typename Sim1::Fp_type>(nqubits);
	sim1.setState(temp);
	sim2.setState(convertState<typename Sim2::Fp_type>(temp));

	std::cout << "Finished constructing benchmark class" << std::endl;
    }

    template<typename Sim1, typename Sim2>
    Results<unsigned, double> Compare<Test::SingleSim, Sim1, Sim2>::impl(
	const std::string & name,
	OneQubitGate<Sim1> fn,
	OneQubitGate<Sim2> gn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;
	
	Results<unsigned, double> res{{"qubit", Sim1::name, Sim2::name}};
	setMeta(name, res);
    
	for (unsigned n = 0; n < nqubits; n++) {
		      
	    qsl::Timer tmr;
		      
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(fn, sim1, n);
	    }
	    tmr.stop();
	    double time1 = tmr.getElapsed();
		    
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(gn, sim2, n);
	    }
	    tmr.stop();
	    double time2 = tmr.getElapsed();

	    printRow(std::to_string(n), time1, time2);
	    res.addRow({n}, {time1, time2});

	}

	return res;
    }

    template<typename Sim1, typename Sim2>
    Results<unsigned, double> Compare<Test::SingleSim, Sim1, Sim2>::impl(
	const std::string & name,
	OneQubitGate<Sim1, typename Sim1::Fp_type> fn,
	OneQubitGate<Sim2, typename Sim2::Fp_type> gn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;

	Results<unsigned, double> res{{"qubit", Sim1::name, Sim2::name}};
	setMeta(name, res);
    
	std::vector<typename Sim1::Fp_type> sim1_phase_list
	    = makeRandomPhases<typename Sim1::Fp_type>(nqubits);
	std::vector<typename Sim2::Fp_type> sim2_phase_list
	    = convertVector<typename Sim2::Fp_type>(sim1_phase_list);
	
	for (unsigned n = 0; n < nqubits; n++) {
		      
	    qsl::Timer tmr;
		      
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(fn, sim1, n, sim1_phase_list[n]);
	    }
	    tmr.stop();
	    double time1 = tmr.getElapsed();
	    
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		std::invoke(gn, sim2, n, sim2_phase_list[n]);
	    }
	    tmr.stop();
	    double time2 = tmr.getElapsed();

	    printRow(std::to_string(n), time1, time2);
	    res.addRow({n}, {time1, time2});
	}

	return res;
    }

    template<typename Sim1, typename Sim2>
    Results<unsigned, double> Compare<Test::SingleSim, Sim1, Sim2>::impl(
	const std::string & name,
	TwoQubitGate<Sim1> fn,
	TwoQubitGate<Sim2> gn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;
    
	Results<unsigned, double> res{{"Control", "Target", Sim1::name, Sim2::name}};
	setMeta(name, res);
    
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned m = n + 1; m < nqubits; m++) {

		qsl::Timer tmr;
	
		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(fn, sim1, n, m);
		}
		tmr.stop();
		double time1 = tmr.getElapsed();
		
		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(gn, sim2, n, m);
		}
		tmr.stop();
		double time2 = tmr.getElapsed();

		std::string col1 = "(" + std::to_string(n) + ", ";
		col1 += std::to_string(m) + ")";
		printRow(col1, time1, time2);
		res.addRow({n, m}, {time1, time2});
	    }
	}

	return res;
    }

    template<typename Sim1, typename Sim2>
    Results<unsigned, double> Compare<Test::SingleSim, Sim1, Sim2>::impl(
	const std::string & name,
	TwoQubitGate<Sim1, typename Sim1::Fp_type> fn,
	TwoQubitGate<Sim2, typename Sim2::Fp_type> gn)
    {
	std::cout << std::endl << "Testing '" << name << "'" << std::endl;

	Results<unsigned, double> res{{"Control", "Target", Sim1::name, Sim2::name}};
	setMeta(name, res);
    
	std::vector<typename Sim1::Fp_type> sim1_phase_list
	    = makeRandomPhases<typename Sim1::Fp_type>(nqubits*nqubits);
	std::vector<typename Sim2::Fp_type> sim2_phase_list
	    = convertVector<typename Sim2::Fp_type>(sim1_phase_list);
    
	for (unsigned n = 0; n < nqubits; n++) {
	    for (unsigned m = n + 1; m < nqubits; m++) {

		qsl::Timer tmr;
	
		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(fn, sim1, n, m, sim1_phase_list[nqubits*n + m]);
		}
		tmr.stop();
		double time1 = tmr.getElapsed();

		tmr.start();
		for (std::size_t i = 0; i < test_len; i++) {
		    std::invoke(gn, sim2, n, m, sim2_phase_list[nqubits*n + m]);
		}
		tmr.stop();
		double time2 = tmr.getElapsed();

		std::string col1 = "(" + std::to_string(n) + ", ";
		col1 += std::to_string(m) + ")";
		printRow(col1, time1, time2);
		res.addRow({n, m}, {time1, time2});
	    }
	}
	return res;
    }

    template<typename Sim1, typename Sim2>
    Results<unsigned, double> Compare<Test::SingleSim, Sim1, Sim2>::measure()
    {
	Results<unsigned, double> res{{"qubit", Sim1::name, Sim2::name}};
	setMeta("measure", res);
    
	// Perform the measurements
	qsl::Timer tmr;

	for (unsigned n = 0; n < nqubits; n++) {

	    tmr.start();
	    for (int i = 0; i < test_len; i++) {
		sim1.measure(n);
	    }
	    tmr.stop();
	    double time1 = tmr.getElapsed();
	
	    tmr.start();
	    for (int i = 0; i < test_len; i++) {
		sim2.measure(n);
	    }
	    tmr.stop();
	    double time2 = tmr.getElapsed();

	    printRow(std::to_string(n), time1, time2);
	    res.addRow({n}, {time1, time2});
	}

	return res;
    }

    template<typename Sim1, typename Sim2>
    Results<std::size_t, double> Compare<Test::SingleSim, Sim1, Sim2>::sampleAll()
    {
	Results<std::size_t, double> r{{"nsamples", Sim1::name, Sim2::name}};
	setMeta("sampleAll", r);
    
	qsl::Timer tmr;

	// Work out the nsamples step size
	std::size_t step = nsamples/20;
    
	for (std::size_t n = 0; n < nsamples; n += step) {
    
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {
		sim1.sampleAll(n);
	    }
	    tmr.stop();
	    double time1 = tmr.getElapsed();
	
	    tmr.start();
	    for (std::size_t i = 0; i < test_len; i++) {	
		sim2.sampleAll(n);
	    }
	    tmr.stop();
	    double time2 = tmr.getElapsed();

	    printRow(std::to_string(n), time1, time2);
	    r.addRow({n}, {time1, time2});
	}
    
	return r;
    }
}
