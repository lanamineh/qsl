#include "qsl.hpp"
#include "meta.hpp"

template <std::floating_point F, bool D, qsl::parallel_t P>
qsl::basic<F, D, P>::basic(unsigned qubits)
    : qubits_{qubits}, size_{1ULL << qubits_}
{
    if constexpr(D) {
	try {
	    state_.resize(size_, 0);
	} catch (const std::bad_alloc &){
	    throw std::runtime_error("Too many qubits to allocate.");
	}
    } else {
	state_.resize(size_, 0);
    }
    
    state_[0] = 1;
}


template <std::floating_point F, bool D, qsl::parallel_t P>
const std::complex<F> & qsl::basic<F, D, P>::operator[](std::size_t index) const
{
    if constexpr(D) {
	if (index >= size_) {
	    throw std::out_of_range("State vector index is out of range.");
	}
    }

    return state_[index];
}


template <std::floating_point F, bool D, qsl::parallel_t P>
qsl::resize<F, D, P>::resize(unsigned qubits)
    : qubits_{qubits}, size_{1ULL << qubits_}
{
    if constexpr(D) {
	try {
	    state_.resize(size_, 0);
	} catch (const std::bad_alloc &){
	    throw std::runtime_error("Too many qubits to allocate.");
	}
    } else {
	state_.resize(size_, 0);
    }
    
    state_[0] = 1;
}


template <std::floating_point F, bool D, qsl::parallel_t P>
const std::complex<F> & qsl::resize<F, D, P>::operator[](std::size_t index) const
{
    if constexpr(D) {
	if (index >= size_) {
	    throw std::out_of_range("State vector index is out of range.");
	}
    }

    return state_[index];
}


template <std::floating_point F, bool D, qsl::parallel_t P>
qsl::number<F, D, P>::number(unsigned qubits)
    : qubits_{qubits}, size_{1ULL << qubits_}
{
    if constexpr(D) {
	try {
	    state_.resize(size_, 0);
	} catch (const std::bad_alloc &){
	    throw std::runtime_error("Too many qubits to allocate.");
	}
    } else {
	state_.resize(size_, 0);
    }
    
    state_[0] = 1;
}


template <std::floating_point F, bool D, qsl::parallel_t P>
const std::complex<F> & qsl::number<F, D, P>::operator[](std::size_t index) const
{
    if constexpr(D) {
	if (index >= size_) {
	    throw std::out_of_range("State vector index is out of range.");
	}
    }

    return state_[index];
}


#include "explicit.tpp"
