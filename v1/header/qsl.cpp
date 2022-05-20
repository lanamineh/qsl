#include "qsl.hpp"
#include "meta.hpp"
#include <sstream>

template <std::floating_point F, bool D>
static inline const std::complex<F> & get_amplitude(const std::vector<std::complex<F>> & state,
						    std::size_t size, std::size_t index)
{
    if constexpr(D) {
	if (index >= size) {
	    std::stringstream ss;
	    ss << "State vector index " << index << " is out of range for state of size " << size
	       << " (in call to operator[]).";
	    throw std::out_of_range(ss.str());
	}
    }

    return state[index];
}


template <std::floating_point F, bool D>
static inline std::vector<std::complex<F>> make_zero_state(std::size_t size)  
{
    if constexpr(D) {
	try {
	    std::vector<std::complex<F>> state(size, 0);
	    state[0] = 1;
	    return state;
	} catch (const std::bad_alloc &){
	    throw std::runtime_error("Too many qubits to allocate.");
	}
    } else {
	std::vector<std::complex<F>> state(size, 0);
	state[0] = 1;
	return state;
    }
}



template <std::floating_point F, bool D, qsl::parallel_t P>
qsl::basic<F, D, P>::basic(unsigned qubits)
    : qubits_{qubits}, size_{1ULL << qubits_}, state_{make_zero_state<F, D>(size_)}
{}


template <std::floating_point F, bool D, qsl::parallel_t P>
const std::complex<F> & qsl::basic<F, D, P>::operator[](std::size_t index) const
{
    return get_amplitude<F, D>(state_, size_, index);
}


template <std::floating_point F, bool D, qsl::parallel_t P>
qsl::resize<F, D, P>::resize(unsigned qubits)
    : qubits_{qubits}, size_{1ULL << qubits_}, state_{make_zero_state<F, D>(size_)}
{}


template <std::floating_point F, bool D, qsl::parallel_t P>
const std::complex<F> & qsl::resize<F, D, P>::operator[](std::size_t index) const
{
    return get_amplitude<F, D>(state_, size_, index);
}


template <std::floating_point F, bool D, qsl::parallel_t P>
qsl::number<F, D, P>::number(unsigned qubits)
    : qubits_{qubits}, size_{1ULL << qubits_}, state_{make_zero_state<F, D>(size_)}
{}


template <std::floating_point F, bool D, qsl::parallel_t P>
const std::complex<F> & qsl::number<F, D, P>::operator[](std::size_t index) const
{
    return get_amplitude<F, D>(state_, size_, index);
}


#include "explicit.tpp"
