#include <vector>
#include <iostream>

/**
 * \brief Base class for quantum simulators
 *
 * Contains only methods which do not depend on any template
 * parameter other than the floating-point precision. The class
 * is augmented with other classes which specify other parts of
 * the behaviour of the simulation.
 *
 */
template<std::floating_point Fp>
class Base
{
protected:
    unsigned num_qubits;
    std::size_t dim;
    std::vector<Fp> state{};
    
public:

    Base(unsigned num_qubits_in)
	: num_qubits{num_qubits_in}, dim{ 1UL << num_qubits }, state(dim)
	{}

    std::size_t dimension() const { return dim; }
    void print() const
	{
	    std::cout << "State on " << num_qubits << " qubits:" << std::endl;
	    for (std::size_t index = 0; index < dim; index++) {
		std::cout << index << ": " << state[index] << std::endl;
	    }
	}  
};

/**
 * \brief Set the state to a particular value
 *
 * This class sets the state sequentially (without OpenMP). A bool is used
 * to specify whether debugging checks should be performed.
 *
 */
template<std::floating_point Fp, bool Debug>
struct SequentialStateSetter : public Base<Fp>
{
    void setState(const std::vector<Fp> & new_state)
    	{
	    if constexpr (Debug == true) {
		if (new_state.size() != Base<Fp>::state.size()) {
		    throw std::logic_error("Cannot set state: wrong dimension");
		}
	    }
	    Base<Fp>::state = new_state;
    	}
};

/**
 * \brief Same set state function but uses OpenMP
 *
 */
template<std::floating_point Fp, bool Debug>
class OmpStateSetter : public Base<Fp>
{
public:
    void setState(const std::vector<Fp> & new_state)
    	{
	    if constexpr (Debug == true) {
		if (new_state.size() != Base<Fp>::state.size()) {
		    throw std::logic_error("Cannot set state: wrong dimension");
		}
	    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for (std::size_t index = 0; index < Base<Fp>::dim; index++) {
		Base<Fp>::state[index] = new_state[index];
	    }
	}
};

template<std::floating_point Fp, bool Debug,
	 template<std::floating_point, bool> class StateSetterPolicy>
class BetterBase : public StateSetterPolicy<Fp, Debug>
{
#ifndef _OPENMP
    static_assert(not std::is_same_v<StateSetterPolicy<Fp>, OmpStateSetter<Fp>>,
		  "You have requested an OpenMP implementation, but _OPENMP "
		  "is not defined. Did you forget -fopenmp?");
#endif
};

template<std::floating_point Fp, bool Debug,
	 template <std::floating_point, bool> class StateSetterPolicy = SequentialStateSetter>
class Generic : public BetterBase<Fp, Debug, StateSetterPolicy>
{
    
};

int main()
{
    Generic<double, false, OmpStateSetter> q{2};
    q.print();

    q.setState({1,0,0});
    q.print();
}

