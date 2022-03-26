#include <vector>
#include <random>
#include <iostream>
#include <complex>

template<typename S>
long long unsigned estimate_size(S sim)
{
    return sizeof(sim) + sim.dim() * sizeof(std::complex<double>);
}

/// This sim does not store the generator
class sim1
{
    unsigned nqubits_;
    std::size_t dim_;
    std::vector<std::complex<double>> amps_;
public:
    sim1(unsigned nqubits)
	: nqubits_{nqubits}, dim_{1UL << nqubits}, amps_(dim_)
	{ }
    std::size_t dim() const { return dim_; }
};

/// This sim stores the generator
class sim2
{
    unsigned nqubits_;
    std::size_t dim_;
    std::vector<std::complex<double>> amps_;
    std::mt19937_64 gen_;
public:
    sim2(unsigned nqubits, long unsigned seed)
	: nqubits_{nqubits}, dim_{1UL << nqubits}, amps_(dim_), gen_{seed}
	{ }  
    std::size_t dim() const { return dim_; }
};


int main()
{
    std::mt19937 g{};
    std::mt19937_64 h{};
    std::cout << "std::mt19937: " << sizeof(g) << ", "
	      << "std::mt19937_64: " << sizeof(h)
	      << std::endl;

    std::cout << "Sizes of simulators" << std::endl;
    for (unsigned n{0}; n < 16; n++) {
	sim1 q1{n};
	sim2 q2{n, 123};
	std::cout << n << " qubits: "
		  << estimate_size(q1) << ", "
		  << estimate_size(q2) << ", "
		  << double(estimate_size(q2))/estimate_size(q1)
		  << std::endl;
    }
}


    
