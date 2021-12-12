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
 * \brief StateSetter base class to provide generic functions
 *
 * This class contains the code for checking the input state. It 
 * is used as a base class for the StateSetter implementations so
 * that they do not have to duplicate the checking code.
 *
 */
template<std::floating_point Fp>
struct StateSetter : public Base<Fp>
{
    void checkInputState(const std::vector<Fp> & new_state) const {

	// Check the state size
	if (new_state.size() != Base<Fp>::state.size()) {
	    throw std::logic_error("Cannot set state: wrong dimension");
	}

	// Check the input state is normalised? etc.
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
struct SequentialStateSetter : public StateSetter<Fp>
{
    void setState(const std::vector<Fp> & new_state)
    	{
	    if constexpr (Debug == true) {
		StateSetter<Fp>::checkInputState(new_state);
	    }
	    Base<Fp>::state = new_state;
    	}
};

/**
 * \brief Same set state function but uses OpenMP
 *
 */
template<std::floating_point Fp, bool Debug>
class OmpStateSetter : public StateSetter<Fp>
{
public:
    void setState(const std::vector<Fp> & new_state)
    	{
	    if constexpr (Debug == true) {
		StateSetter<Fp>::checkInputState(new_state);		
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
    static_assert(not std::is_same_v<StateSetterPolicy<Fp,Debug>,
		  OmpStateSetter<Fp,Debug>>,
		  "You have requested an OpenMP implementation, but _OPENMP "
		  "is not defined. Did you forget -fopenmp?");
#endif
};

/**
 * \brief A collection of all the settings required for sequential operation
 *
 */
struct Sequential
{
    template<std::floating_point Fp, bool Debug>
    using StateSetterPolicy = SequentialStateSetter<Fp, Debug>;
};

struct Omp
{
    template<std::floating_point Fp, bool Debug>
    using StateSetterPolicy = OmpStateSetter<Fp, Debug>;
};

template<std::floating_point Fp, bool Debug = false, typename SeqPar = Sequential>
class Generic : public BetterBase<Fp, Debug, SeqPar::template StateSetterPolicy>
{
    
};


int main()
{
    Generic<double, true, Sequential> q{2};
    q.print();

    q.setState({1,0,0,0});
    q.print();
}

