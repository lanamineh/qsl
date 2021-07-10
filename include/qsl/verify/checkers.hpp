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
 * \brief Checkers for use with the Verify class
 *
 */

#ifndef CHECKERS_HPP
#define CHECKERS_HPP

#include <memory>
#include <vector>
#include <cstddef>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>

namespace qsl {

    /**
     * \brief Simulator checker object concept
     *
     * All simulator checkers should conform to this specification.
     * A simulator checker must
     *
     *   1) be templated on typename two qsl::Simulator types, Sim1 and
     *      Sim2, the ones that are being checked
     *   2) be default constructible
     *   3) contain a bind(std::unique_ptr<Sim1>&, std::unique_ptr<Sim2>&)
     *      function to attach two simulators to the checker. The simulators
     *      are assumed to be initialised in equal states according to the
     *      StateGen class
     *   4) contain a checkAll() member function, which performs checks on
     *      the simulators.
     *
     * The checker can optionally contain a configureChecker method with
     * any parameter list
     * 
     */
    template<typename T, typename Sim1, typename Sim2>
    concept SimChecker = std::is_default_constructible<T>::value
	&& qsl::Simulator<Sim1>
	&& qsl::Simulator<Sim2> 
    	&& requires(T t, std::unique_ptr<Sim1> p1, std::unique_ptr<Sim2> p2)
    {
    	// Required member functions
    	t.bind(p1,p2);
    };
    
    /**
     * \brief Function to calculate 2-sided z-value
     */
    double getZValue(double ci)
    {
	// Calculate the confidence interval scale
	double p = 0.5 + ci/2;

	// This objective function will be zero if \Phi(x) = p
	auto obj = [=](double x) {
	    return 0.5 * std::erfc(-x * M_SQRT1_2) - p;
	};

	// Secant method
	double x0 = 0, x1 = 1, x2 = 0, f0 = 0, f1 = 1;
	f0 = obj(x0);
    
	for (int i = 0; i < 20; i++) {
	    if (std::abs(f1) < 1e-5) {
		break;
	    }
	    f1 = obj(x1);	
	
	    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
	    x0 = x1;
	    f0 = f1;
	    x1 = x2;
	}

	return x1;
    }

    /**
     * \brief Compute a binomial confidence interval
     * 
     */
    class BinomialCI
    {
	double halfWidth(double p, double level, unsigned nsamples) const
	    {
		double z_value = getZValue(level);
		return z_value * std::sqrt((p * (1-p))/nsamples);
	    }
    
    public:
	const double center;
	const double width;
	const double bottom;
	const double top;
	const double level;
	BinomialCI(double p, double level_in, std::size_t nsamples)
	    : center{p}, width{2*halfWidth(p, level_in, nsamples)},
	      bottom{center-width/2}, top{center+width/2},
	      level{level_in}
	    {}
    };

    /// Return true if the confidence intervals overlap
    bool overlap(BinomialCI a, BinomialCI b)
    {
	return std::abs(a.center - b.center) < 0.5*(a.width + b.width);
    }

    /// Compute the maximum relative error given two confidence intervals
    double maxRelativeError(BinomialCI a, BinomialCI b)
    {
	if (not overlap(a,b)) {
	    throw std::logic_error("Confidence intervals do not overlap");
	}
	// Compute largest possible error given confidence intervals
	double max_p = std::max(a.top, b.top);
	double min_p = std::min(a.bottom, b.bottom);
	///\todo What if both are zero?
	double max_relative_error = (max_p - min_p)/max_p;		
	return max_relative_error;
    
    };

    /// Copy and paste this class to make a Checker
    template<typename Sim1, typename Sim2>
    class ExampleChecker
    {
	std::unique_ptr<Sim1> sim1;
	std::unique_ptr<Sim2> sim2;

    public:

	/// Do the check
	void checkAll()
	    {
		if (not sim1 or not sim2) {
		    throw std::logic_error(
			"Cannot run check before binding simulators");
		}

		std::cout << "Put the check here" << std::endl;
	    }

	/// Set up the checker here
	void configureChecker(int val)
	    {
		std::cout << "Doing something different: " << val << std::endl;
	    }

    
	/// Attach simulator objects to this checker
	void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	    {
		std::swap(sim1, p1);
		std::swap(sim2, p2);
	    }
    };

    /**
     * \brief Check the prob function 
     *
     * The prob function returns the probability of getting a particular
     * outcome from a measurement. The prob function takes a qubit index, 
     * and an outcome, and returns the probability of getting that outcome
     * when measuring the qubit. The internal state of the simulator is
     * not modified by the function.
     *
     * The prob function is deterministic, so checking it amounts to checking
     * that the results from the two simulators are equal.
     *
     */
    template<qsl::Simulator Sim1, qsl::Simulator Sim2>
    class ProbChecker
    {
	std::unique_ptr<Sim1> sim1;
	std::unique_ptr<Sim2> sim2;

    public:

	/// Do the check
	void checkAll()
	    {
		if (not sim1 or not sim2) {
		    throw std::logic_error(
			"Cannot run check before binding simulators");
		}

		const unsigned nqubits = sim1->getNumQubits();
		for (unsigned n = 0; n < nqubits; n++) {
		    double p1_out0 = sim1->prob(n, 0);
		    double p1_out1 = sim1->prob(n, 1);
		    double p2_out0 = sim2->prob(n, 0);
		    double p2_out1 = sim2->prob(n, 1);

		    std::cout << "n = " << n << ": "
			      << "prob 0 diff = " << p1_out0 - p2_out0
			      << ", prob 1 diff = " << p1_out1 - p2_out1
			      << std::endl;
		}
	    }

	/// Set up the checker here
	void configureChecker(int val)
	    {
	    }

    
	/// Attach simulator objects to this checker
	void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	    {
		std::swap(sim1, p1);
		std::swap(sim2, p2);
	    }
    };

    /**
     * \brief Check the measure function
     *
     * This class checks the measure function. The measure function takes
     * one argument, which is the qubit to be measured, and return either
     * zero or one which is the measurement outcome. In addition, the
     * state is collapsed according to the measurement outcome.
     * 
     */
    template<qsl::Simulator Sim1, qsl::Simulator Sim2>
    class MeasureChecker
    {
	std::size_t nsamples;
	double ci; // confidence level
	std::unique_ptr<Sim1> sim1;
	std::unique_ptr<Sim2> sim2;

    public:

	MeasureChecker() : nsamples(100), ci{0.95} {}
    
	/**
	 * \brief Check that the measure functions are the same
	 *
	 *
	 */
	void checkAll()
	    {
		if (not sim1 or not sim2) {
		    throw std::logic_error(
			"Cannot run check before binding simulators");
		}

		const unsigned nqubits = sim1->getNumQubits();

		// Store copies of the simulators to repeat the test
		Sim1 sim1_copy{ sim1->getNumQubits() };
		sim1_copy.setState(sim1->getState());
		Sim2 sim2_copy{ sim2->getNumQubits() };
		sim2_copy.setState(sim2->getState());
		
		// Repeat the test 
		
		for (unsigned n = 0; n < nqubits; n++) {

		    // Vectors to store the repeated measured results
		    std::vector<int> outcomes1;
		    std::vector<int> outcomes2;
		    
		    // Repeat the test nsamples times
		    for (std::size_t k = 0; k < nsamples; k++) {

			// Reset the simulators back to the copy
			sim1->setState(sim1_copy.getState());
			sim2->setState(sim2_copy.getState());
			
			// Call the measure function and store the result
			outcomes1.push_back(sim1->measure(n));
			outcomes2.push_back(sim2->measure(n));

			// If the outcomes are the same, check that
			// the state vectors are also the same
			if (outcomes1.back() == outcomes2.back()) {
			    auto state1 = sim1->getState();
			    auto state2 = sim2->getState();
			    double distance = fubiniStudy(state1, state2);
			    std::cout << "Distance = " << distance << std::endl; 
			    ///\todo Count up how many times the distances
			    /// are zero
			}
			
		    }

		    // Find the total number of ones by adding up all the
		    // outcomes (the zeros contribute nothing). Then
		    // compute the expected value of 1 (equal to Binomial p).
		    // 
		    unsigned sum1 = std::accumulate(std::begin(outcomes1),
						    std::end(outcomes1),
						    0);
		    unsigned sum2 = std::accumulate(std::begin(outcomes2),
						    std::end(outcomes2),
						    0);
		    double p1 = static_cast<double>(sum1)/nsamples;
		    double p2 = static_cast<double>(sum2)/nsamples;
		
		    // Compute (half) the confidence interval width
		    BinomialCI ci1{p1, ci, nsamples};
		    BinomialCI ci2{p2, ci, nsamples};
	       		
		    // Check if confidence intervals overlap
		    std::cout << "n = " << n << ": "
			      << "p1 = " << p1 << ", p2 = " << p2 << ", "; 
		    if (overlap(ci1, ci2)) {
			// The ranges overlap
			std::cout << "error < "
				  << 100*maxRelativeError(ci1,ci2) << "%"
				  << " (" << 100*ci*ci << "% level)"
				  << std::endl; 
		    } else {
			// Ranges do not overlap. If p1 is in range 1
			// with prob ci, and p2 is in range 2 with
			// prob ci, then p1 != p2 with prob ci^2
			std::cout << "p1 != p2 "
				  << " (" << 100*ci*ci << "% level)"
				  << std::endl;
		    }
		}
		
	    }

	/// Set up the checker here
	///\todo Work out the right name for the variable ci (confidence level?)
	void configureChecker(std::size_t nsamples_in, double ci_in)
	    {
		std::cout << "Setting sample number (nsamples) = " << nsamples_in
			  << std::endl;
		nsamples = nsamples_in;

		ci = ci_in;
		std::cout << "Setting confidence limit to = "
			  << ci*100 << "% (using z = " << getZValue(ci) << ")"
			  << std::endl;
	    }

    
	/// Attach simulator objects to this checker
	void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	    {
		std::swap(sim1, p1);
		std::swap(sim2, p2);
	    }
    };

    
    /// Check the sample function
    template<qsl::Simulator Sim1, qsl::Simulator Sim2>
    class SampleChecker
    {
	std::size_t nsamples;
	double ci; // confidence level
	std::unique_ptr<Sim1> sim1;
	std::unique_ptr<Sim2> sim2;

    public:

	SampleChecker() : nsamples(100), ci{0.95} {}
    
	/**
	 * \brief Check that the sample functions are the same
	 *
	 * The functions calls sample of the two simulators with
	 * large sample sizes, estimates the probability of getting
	 * a 1 for each target qubit, then calculates a confidence
	 * interval about each estimate. If the confidence intervals
	 * overlap the the test is a success, and the error is printed.
	 * If the intervals do not overlap, then the test fails.
	 *
	 * The confidence interval is calculated using the formula for
	 * confidence interval on the binomial distribution confidence
	 * intervals wikipedia page.
	 *
	 * \todo For this to be legitimate with low sample numbers, it
	 * is necessary to repeatedly sample the sample function (which 
	 * is not what happens at the moment).
	 *
	 */
	void checkAll()
	    {
		if (not sim1 or not sim2) {
		    throw std::logic_error(
			"Cannot run check before binding simulators");
		}

		const unsigned nqubits = sim1->getNumQubits();

		for (unsigned n = 0; n < nqubits; n++) {
		    // Get the samples
		    std::vector<std::size_t> sample1 = sim1->sample(n, nsamples);
		    std::vector<std::size_t> sample2 = sim2->sample(n, nsamples);

		    // Compute the expected value of 1 (equal to Binomial p)
		    double p1 = static_cast<double>(sample1[1])/nsamples;
		    double p2 = static_cast<double>(sample2[1])/nsamples;
		
		    // Compute (half) the confidence interval width
		    BinomialCI ci1{p1, ci, nsamples};
		    BinomialCI ci2{p2, ci, nsamples};
	       		
		    // Check if confidence intervals overlap
		    std::cout << "n = " << n << ": "
			      << "p1 = " << p1 << ", p2 = " << p2 << ", "; 
		    if (overlap(ci1, ci2)) {
			// The ranges overlap
			std::cout << "error < "
				  << 100*maxRelativeError(ci1,ci2) << "%"
				  << " (" << 100*ci*ci << "% level)"
				  << std::endl; 
		    } else {
			// Ranges do not overlap. If p1 is in range 1
			// with prob ci, and p2 is in range 2 with
			// prob ci, then p1 != p2 with prob ci^2
			std::cout << "p1 != p2 "
				  << " (" << 100*ci*ci << "% level)"
				  << std::endl;
		    }

		}
	    
	    }

	/// Set up the checker here
	///\todo Work out the right name for the variable ci (confidence level?)
	void configureChecker(std::size_t nsamples_in, double ci_in)
	    {
		std::cout << "Setting sample number (nsamples) = " << nsamples_in
			  << std::endl;
		nsamples = nsamples_in;

		ci = ci_in;
		std::cout << "Setting confidence limit to = "
			  << ci*100 << "% (using z = " << getZValue(ci) << ")"
			  << std::endl;
	    }

    
	/// Attach simulator objects to this checker
	void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	    {
		std::swap(sim1, p1);
		std::swap(sim2, p2);
	    }
    };

    /**
     * \brief Check the sampelAll function
     *
     * The sampleAll function simulates repeated measurements of all the qubits,
     * but without actually collapsing and re-preparing the state vector. The
     * results are returned in a standard map, which associates each outcome
     * that occured with its frequency in the sample.
     *
     * The checker works but taking the sampled results from two simulators and
     * comparing the relative frequency of each result. Each outcome is treated
     * separately as binomially distributed (the result was either the outcome
     * or it wasn't), which allows confidence intervals to be calculated for the
     * relative frequency of each outcome.
     *
     */
    template<qsl::Simulator Sim1, qsl::Simulator Sim2>
    class SampleAllChecker
    {
	std::size_t nsamples;
	double ci; // confidence level
	std::unique_ptr<Sim1> sim1;
	std::unique_ptr<Sim2> sim2;

    public:

	SampleAllChecker() : nsamples(100), ci{0.95} {}
    
	/**
	 * \brief Check that the sample functions are the same
	 *
	 * The function prints the average relative error over all overlapping
	 * confidence intervals for relative frequencies. It also prints the
	 * proportion of relative frequencies whose confidence intervals overlap,
	 * and the proportion of outcomes which occur in both sample1 and sample2.
	 *
	 * \todo Need to make the proportion of common outcomes symmetrical with
	 * respect to sample2.
	 */
	void checkAll()
	    {
		if (not sim1 or not sim2) {
		    throw std::logic_error(
			"Cannot run check before binding simulators");
		}

		// Get the samples
		std::map<std::size_t, std::size_t> sample1
		    = sim1->sampleAll(nsamples);
		std::map<std::size_t, std::size_t> sample2
		    = sim2->sampleAll(nsamples);

		// Relative error for outcomes whose confidence
		// intervals overlap
		std::map<std::size_t, double> relative_errors;

		// Number of outcomes where the confidence intervals
		// did not overlaps
		unsigned num_non_overlaps = 0;

		// Number of outcomes in both sample1 and sample2
		unsigned common_outcomes = 0;

		// Calculate the relative frequencies
		for(const auto & [outcome,frequency1] : sample1) {

		    // Compute confidence interval for outcome in sample1
		    double p1 = static_cast<double>(frequency1)/nsamples;
		    BinomialCI ci1{p1, ci, nsamples};

		    try {
			// Compute confidence interval for outcome in sample2
			std::size_t frequency2 = sample2.at(outcome);
			double p2 = static_cast<double>(frequency2)/nsamples;
			BinomialCI ci2{p2, ci, nsamples};

			// If the previous part worked, then outcome is
			// common to both sample1 and sample2
			common_outcomes++;
			if (overlap(ci1, ci2)) {
			    // The ranges overlap
			    relative_errors[outcome] = maxRelativeError(ci1,ci2);
			} else {
			    // Ranges do not overlap. If p1 is in range 1
			    // with prob ci, and p2 is in range 2 with
			    // prob ci, then p1 != p2 with prob ci^2
			    num_non_overlaps++;
			}
			

		    } catch (const std::out_of_range & e) {
			//std::cout << "Outcome not present in sample2"
			//	  << std::endl;
		    }
		}

		// Calculate the mean relative error
		double mean_relative_error = 0;
		for (const auto item : relative_errors) {
		    mean_relative_error += item.second;
		}
		mean_relative_error /= relative_errors.size();

		// Compute overlap rate, the proportion of outcomes
		// with confidence intervals that overlapped in
		// sample1 and sample2
		double overlap_rate
		    = 1 - static_cast<double>(num_non_overlaps)/common_outcomes;

		// Compute the proportion of outcomes common to both
		// sample1 and sample2.
		///\todo Fix so it is symmetrical with respect to sample2
		double common_rate
		    = static_cast<double>(common_outcomes)/sample2.size();
		
		std::cout << "Mean error = "
			  << 100*mean_relative_error << "%"
			  << ", overlap = " << 100*overlap_rate << "%"
			  << ", common = " << 100*common_rate << "%"
			  << std::endl;
		
	    }

	/// Set up the checker here
	///\todo Work out the right name for the variable ci (confidence level?)
	void configureChecker(std::size_t nsamples_in, double ci_in)
	    {
		std::cout << "Setting sample number (nsamples) = " << nsamples_in
			  << std::endl;
		nsamples = nsamples_in;

		ci = ci_in;
		std::cout << "Setting confidence limit to = "
			  << ci*100 << "% (using z = " << getZValue(ci) << ")"
			  << std::endl;
	    }

    
	/// Attach simulator objects to this checker
	void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	    {
		std::swap(sim1, p1);
		std::swap(sim2, p2);
	    }
    };


    /**
     * \brief Check gates
     *
     * This class contains the routines that are used to check 
     * gates. There is no checkAll method, so it cannot be used
     * as an argument in the Verify class. However, classes 
     * derived from this one that implement checkAll can be used
     * in the Verify class.
     *
     * The class does nothing other than implement gates. It does
     * not modify the simulators in any other way (e.g. setting 
     * states, etc.)
     *
     * \todo Find a way to control the phases, etc.
     */
    template<qsl::Simulator Sim1, qsl::Simulator Sim2>
    class GateChecker
    {
	std::unique_ptr<Sim1> sim1;
	std::unique_ptr<Sim2> sim2;

    public:
    
	/// Implementation for one qubit gates with no arguments
	void check(OneQubitGate<Sim1> fn, OneQubitGate<Sim2> gn)
	    {
		unsigned nqubits = sim1->getNumQubits();
	    
		// Loop through every pair of sim1
		for (unsigned n = 0; n < nqubits; n++) {

		    // equate the states
		    sim2->setState(sim1->getState());

		    std::cout << "Qubit no " << n << ": ";
		
		    // Apply the gates stored in member pointers fn and gn
		    std::invoke(fn, sim1, n);
		    std::invoke(gn, sim2, n);
		
		    ///\todo Should be float?
		    double distance = fubiniStudy(sim1->getState(),
						  sim2->getState());

		    std::cout << "Distance = " << distance << std::endl;
		}
	    }

	/// Checkementation for single qubit gates with one double argument
	void check(OneQubitGate<Sim1, typename Sim1::Fp_type> fn,
		   OneQubitGate<Sim2, typename Sim2::Fp_type> gn)
	    {
		unsigned nqubits = sim1->getNumQubits();
		std::vector<typename Sim1::Fp_type> sim1_phase_list
		    = makeRandomPhases<typename Sim1::Fp_type>(nqubits);
		std::vector<typename Sim2::Fp_type> sim2_phase_list
		    = convertVector<typename Sim2::Fp_type>(sim1_phase_list);
    
		// Loop through every pair of sim1
		for (unsigned n = 0; n < nqubits; n++) {

		    // equate the states
		    sim2->setState(sim1->getState());

		    std::cout << "Qubit no " << n << ": ";
		    // Apply the gates stored in member pointers fn and gn
		    std::invoke(fn, sim1, n, sim1_phase_list[n]);
		    std::invoke(gn, sim2, n, sim2_phase_list[n]);
		
		    ///\todo Should be float?
		    double distance = fubiniStudy(sim1->getState(),
						  sim2->getState());
		    std::cout << "Distance = " << distance << std::endl;
		}
	    }

	/// Checkementation for two qubit gates with no arguments
	void check(TwoQubitGate<Sim1> fn, TwoQubitGate<Sim2> gn)
	    {
		unsigned nqubits = sim1->getNumQubits();
    
		// Loop through every pair of sim1
		for (unsigned n = 0; n < nqubits; n++) {
		    for (unsigned m = 0; m < nqubits; m++) {

			// equate the states
			sim2->setState(sim1->getState());

			if (n != m) {
			    std::cout << "Sim1 " << n  << " and " << m << ": ";

			    // Apply the gates stored in member pointers fn and gn
			    std::invoke(fn, sim1, n, m);
			    std::invoke(gn, sim2, n, m);
		
			    double distance = fubiniStudy(sim1->getState(),
							  sim2->getState());
			    std::cout << "Distance = " << distance << std::endl;
			}
		    }
		}
	    }

	/// Checkementation for two qubit gates with one double argument
	void check(TwoQubitGate<Sim1, typename Sim1::Fp_type> fn,
		   TwoQubitGate<Sim2, typename Sim2::Fp_type> gn)
	    {
	
		unsigned nqubits = sim1->getNumQubits();
		//std::vector<double> phase_list = makeRandomPhases(nqubits*nqubits);
		std::vector<typename Sim1::Fp_type> sim1_phase_list
		    = makeRandomPhases<typename Sim1::Fp_type>(nqubits*nqubits);
		std::vector<typename Sim2::Fp_type> sim2_phase_list
		    = convertVector<typename Sim2::Fp_type>(sim1_phase_list);
    
		// Loop through every pair of sim1
		for (unsigned n = 0; n < nqubits; n++) {
		    for (unsigned m = 0; m < nqubits; m++) {

			// equate the states
			sim2->setState(sim1->getState());

			if (n != m) {
			    std::cout << "Sim1 " << n  << " and " << m << ": ";

			    // Apply the gates stored in member pointers fn and gn
			    std::invoke(fn, sim1, n, m, sim1_phase_list[nqubits*n + m]);
			    std::invoke(gn, sim2, n, m, sim2_phase_list[nqubits*n + m]);

			    ///\todo Should be float or Fp_type?
			    double distance = fubiniStudy(sim1->getState(),
							  sim2->getState());
			    std::cout << "Distance = " << distance << std::endl;
			}
		    }
		}
	    }


	/**
	 * \brief Attach simulator objects to this checker
	 */
	void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	    {
		std::swap(sim1, p1);
		std::swap(sim2, p2);
	    }

	void configureChecker(int val)
	    {
		std::cout << "Doing something different: " << val << std::endl;
	    }
    
    }; 

    /// Example gate checker (check all gates)
    template<HasAllGates Sim1, HasAllGates Sim2>
    class DefaultGateChecker : public GateChecker<Sim1,Sim2>
    {
    public:
	void checkAll() {

	    std::cout << "Checking pauliX" << std::endl;
	    this->check(&Sim1::pauliX, &Sim2::pauliX);
	    std::cout << "Checking rotateX" << std::endl;
	    this->check(&Sim1::rotateX, &Sim2::rotateX);
	    std::cout << "Checking phase" << std::endl;
	    this->check(&Sim1::phase, &Sim2::phase);
	    std::cout << "Checking controlNot" << std::endl;
	    this->check(&Sim1::controlNot, &Sim2::controlNot);
	    std::cout << "Checking controlPhase" << std::endl;
	    this->check(&Sim1::controlPhase, &Sim2::controlPhase);
	    std::cout << "Checking swap" << std::endl;
	    this->check(&Sim1::swap, &Sim2::swap);
	    std::cout << "Checking hadamard" << std::endl;
	    this->check(&Sim1::hadamard, &Sim2::hadamard);
	    std::cout << "Checking controlZ" << std::endl;
	    this->check(&Sim1::controlZ, &Sim2::controlZ);
	    
	}
    };

    /**
     * \brief Check number preserved gates
     *
     * Pass this class to the Verify class to check number
     * preserving gates.
     *
     */
    template<HasNPGates Sim1, HasNPGates Sim2>
    class NPGateChecker : public GateChecker<Sim1,Sim2>
    {
    public:
	void checkAll() {
	    std::cout << "Checking phase" << std::endl;
	    this->check(&Sim1::phase, &Sim2::phase);
	    std::cout << "Checking controlPhase" << std::endl;
	    this->check(&Sim1::controlPhase, &Sim2::controlPhase);
	    std::cout << "Checking swap" << std::endl;
	    this->check(&Sim1::swap, &Sim2::swap);
	    std::cout << "Checking controlZ" << std::endl;
	    this->check(&Sim1::controlZ, &Sim2::controlZ);
	    
	}
    };

    /**
     * \brief Check measurement related functions
     *
     * The purpose of this class is to contain checks for all the
     * measurement functions. It is not supposed to be used directly
     * (it doesn't have a checkAll function). Instead, classes
     * derived from this class can call a subset of the measurement 
     * checks depending on what is required.
     *
     */
    template<qsl::Simulator Sim1, qsl::Simulator Sim2>
    class PostselectChecker
    {
	std::unique_ptr<Sim1> sim1;
	std::unique_ptr<Sim2> sim2;

    public:

	/**
	 * \brief Check postselection
	 *
	 * For concepts: requires copy assignment, copy, postselect and getState.
	 */
	void checkAll()
	    {
		if (not sim1 or not sim2) {
		    throw std::logic_error(
			"Cannot run check before binding simulators");
		}
	    
		unsigned nqubits = sim1->getNumQubits();
		std::cout << "============================================"
			  << std::endl;
		std::cout << "Verifying postselect on " << nqubits << " qubits."
			  << std::endl;

		std::cout << "Printing probabilities of collapsing to outcome"
			  << std::endl;
    
		for (unsigned n = 0; n < nqubits; n++) {

		    // Make local copy
		    Sim1 q1(*sim1); // Make a local copy
		    Sim2 q2(*sim2); // Make a local copy
	
		    // Postselect on 0
		    double qubits_0 = q1.postselect(n, 0);
		    double quest_0 = q2.postselect(n, 0);
		    double norm_0 = norm(q1.getState());
		    double distance_0 = fubiniStudy(q1.getState(), q2.getState());

		    // Reset the state vectors
		    q1 = *sim1;
		    q2 = *sim2;

		    // Postselect on 1
		    double qubits_1 = q1.postselect(n, 1);
		    double quest_1 = q2.postselect(n, 1);
		    double norm_1 = norm(q1.getState());
		    double distance_1 = fubiniStudy(q1.getState(), q2.getState());
	
		    std::cout << "Qubit " << n
			      << ": Qubits = [" << qubits_0 << ", " << qubits_1 << "]"
			      << ", Quest = [" << quest_0 << ", " << quest_1 << "]"
			      << std::endl
			      << "Qubits = [" << distance_0 << ", " << distance_1 << "]"
			      << ", Qubits norms = [" << norm_0 << ", "
			      << norm_1 << "]" << std::endl;
		}
	    }

	/**
	 * \brief Attach simulator objects to this checker
	 */
	void bind(std::unique_ptr<Sim1> & p1, std::unique_ptr<Sim2> & p2)
	    {
		std::swap(sim1, p1);
		std::swap(sim2, p2);
	    }
    
    
    };

}
    
#endif
