#include <vector>
#include <concepts>
#include <stdexcept>

template <std::floating_point F>
class Sim
{
public:
    unsigned nqubits;
    std::size_t dim;
    std::vector<F> state;
    
    Sim<F>(unsigned nq) : nqubits{nq}, dim{1ULL << nq}
	{
	    if (nqubits > 10) {
		throw std::logic_error("Too many qubits to simulate");
	    }
	    else {
		state.resize(dim, 0);
	    }
	    
	};
};


template <std::floating_point F, bool D>
class Sim2
{
public:
    constexpr bool debug() { return D; }
    
    unsigned nqubits;
    std::size_t dim;
    std::vector<F> state;
    
    Sim2<F, D>(unsigned nq) : nqubits{nq}, dim{1ULL << nq}
	{
	    if (nqubits > 10) {
		throw std::logic_error("Too many qubits to simulate");
	    }
	    else {
		state.resize(dim, 0);
	    }
	    
	};
};

