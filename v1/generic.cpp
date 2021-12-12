#include <vector>
#include <iostream>

// struct SequentialSetStatePolicy
// {
//     template<std::floating_point Fp>
//     static void setState([[maybe_unused]] std::size_t dim,
// 			 std::vector<Fp> & state,
// 			 const std::vector<Fp> & new_state)
// 	{
// 	    state = new_state;
// 	}
// };

// struct OmpSetStatePolicy
// {
//     template<std::floating_point Fp>
//     static void setState(std::size_t dim,
// 			 std::vector<Fp> & state,
// 			 const std::vector<Fp> & new_state)
// 	{
// #pragma omp parallel for
// 	    for (std::size_t index= 0; index < dim; index++) {
// 		state[index] = new_state[index];
// 	    }
// 	}
// };



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

    // void setState(const std::vector<Fp> & new_state)
    // 	{
    // 	    SetStatePolicy::template setState<Fp>(dim, state, new_state);
    // 	}
    
};

template<std::floating_point Fp, bool Debug = false>
struct SequentialStateSetter : public Base<Fp>
{
    void setState(const std::vector<Fp> & new_state)
    	{
	    Base<Fp>::state = new_state;
    	}
};

template<std::floating_point Fp>
struct SequentialStateSetter<Fp, true> : public Base<Fp>
{
    void setState(const std::vector<Fp> & new_state)
    	{
	    if (new_state.size() != Base<Fp>::state.size()) {
		throw std::logic_error("Cannot set state: wrong dimension");
	    }
	    Base<Fp>::state = new_state;
    	}    
};

template<std::floating_point Fp, bool Debug>
class OmpStateSetter : public Base<Fp>
{
public:
    void setState(const std::vector<Fp> & new_state)
    	{
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for (std::size_t index = 0; index < Base<Fp>::dim; index++) {
		Base<Fp>::state[index] = new_state[index];
	    }
	}
};

template<std::floating_point Fp>
class OmpStateSetter<Fp, true> : public Base<Fp>
{
public:
    void setState(const std::vector<Fp> & new_state)
    	{
	    if (new_state.size() != Base<Fp>::state.size()) {
		throw std::logic_error("Cannot set state: wrong dimension");
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
    Generic<double, true, OmpStateSetter> q{2};
    q.print();

    q.setState({1,0,0});
    q.print();
}

