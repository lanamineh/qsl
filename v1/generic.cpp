#include <vector>
#include <iostream>

struct SequentialSetStatePolicy
{
    template<std::floating_point Fp>
    static void setState([[maybe_unused]] std::size_t dim,
			 std::vector<Fp> & state,
			 const std::vector<Fp> & new_state)
	{
	    state = new_state;
	}
};

struct OmpSetStatePolicy
{
    template<std::floating_point Fp>
    static void setState(std::size_t dim,
			 std::vector<Fp> & state,
			 const std::vector<Fp> & new_state)
	{
#pragma omp parallel for
	    for (std::size_t index= 0; index < dim; index++) {
		state[index] = new_state[index];
	    }
	}
};



template<std::floating_point Fp, class SetStatePolicy>
class Base
{
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

    void setState(const std::vector<Fp> & new_state)
	{
	    SetStatePolicy::template setState<Fp>(dim, state, new_state);
	}
    
};


// struct Omp_Yes
// {
//     using SetStatePolicy = OmpSetStatePolicy;
// };

struct Omp_No
{
    using SetStatePolicy = SequentialSetStatePolicy;
};


template<std::floating_point Fp, class Omp>
class Generic : public Base<Fp, typename Omp::SetStatePolicy>
{
    
};

int main()
{
    Generic<double, Omp_No> q{2};
    q.print();

    q.setState({1,0,0,1});
    q.print();
}

