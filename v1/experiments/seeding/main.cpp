#include <vector>
#include <random>
#include <iostream>
#include <complex>
#include <fstream>

/**
 * \brief Get a random seed
 *
 * According to the page below, it is best to avoid std::random_device,
 * because it is not guaranteed to be random:
 * 
 * https://stackoverflow.com/questions/45069219/
 *   how-to-succinctly-portably-and-thoroughly-
 *   seed-the-mt19937-prng
 *
 * This method uses /dev/urandom to get a random seed.
 * 
 */
template<std::integral T>
static T sysrandom()
{
    T seed;
    char * buffer = reinterpret_cast<char*>(&seed);
    std::ifstream stream("/dev/urandom",
			 std::ios_base::binary | std::ios_base::in);
    stream.read(buffer, sizeof(seed));
    return seed;
}

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
    std::complex<double> get() const { return amps_[0]; }
};

/// This sim stores the generator
class sim2
{
    unsigned nqubits_;
    std::size_t dim_;
    std::vector<std::complex<double>> amps_;
    std::mt19937_64 gen_;
public:
    sim2(unsigned nqubits, long unsigned seed = sysrandom<std::uint_fast64_t>())
	: nqubits_{nqubits}, dim_{1UL << nqubits}, amps_(dim_), gen_{seed}
	{ }  
    std::size_t dim() const { return dim_; }
    std::complex<double> get() const { return amps_[0]; }
};


int main()
{
    // std::mt19937 g{};
    // std::mt19937_64 h{};
    // std::cout << "std::mt19937: " << sizeof(g) << ", "
    // 	      << "std::mt19937_64: " << sizeof(h)
    // 	      << std::endl;

    // std::cout << "Sizes of simulators" << std::endl;
    // for (unsigned n{0}; n < 16; n++) {
    // 	sim1 q1{n};
    // 	sim2 q2{n, 123};
    // 	std::cout << n << " qubits: "
    // 		  << sizeof(q1) << ", "
    // 		  << estimate_size(q2) << ", "
    // 		  << double(estimate_size(q2))/estimate_size(q1)
    // 		  << std::endl;
    // }

    std::complex<double> res{0,0};
    for (std::size_t n{0}; n < 10000000; n++) {
    std::vector<double> abc(1100);
	sim2 s{1};
        res += s.get();
	sim2 s2{s};
	res += (s.get() + s2.get());
    }
    std::cout << res << std::endl; 
}


    
