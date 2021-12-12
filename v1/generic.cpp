#include <vector>
#include <iostream>

struct Omp_Yes{};
struct Omp_No{};

template<std::floating_point Fp, typename Omp>
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
	    static_assert(std::is_same_v<Omp, Omp_No>,
			  "You are trying to use an OpenMP implementation, "
			  "but _OPENMP is not defined. Did you forget "
			  "-fopenmp?");
	    state = new_state;
	}
};

#ifdef _OPENMP
template<std::floating_point Fp>
void Base<Fp, Omp_Yes>::setState(const std::vector<Fp> & new_state) {
#pragma omp parallel for
    for (std::size_t index = 0; index < dim; index++) {
	state[index] = new_state[index];
    }
}
#endif

template<std::floating_point Fp>
void Base<Fp, Omp_Yes>::setState(const std::vector<Fp> & new_state) {
#pragma omp parallel for
    for (std::size_t index = 0; index < dim; index++) {
	state[index] = new_state[index];
    }
}

template<std::floating_point Fp, typename Omp>
class Generic : public Base<Fp, Omp>
{
    
};

int main()
{
    Generic<double, Omp_Yes> q{2};
    q.print();

    q.setState({1,0,0,1});
    q.print();
}

