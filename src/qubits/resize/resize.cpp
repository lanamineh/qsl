/*
 *  Authors: Lana Mineh and John Scott
 *  Copyright 2021 Phasecraft Ltd. and John Scott
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */
 
/**
 * \file resize.cpp
 * \brief Contains the implementation of constructors and 
 *        set/get methods of the Resize qubits class.
 */


#include "qsl/qubits.hpp"
#include "qsl/utils/misc.hpp"
#include "qsl/utils/quantum.hpp"
#include <iostream>

template<>
const std::string qsl::Qubits<qsl::Type::Resize, double>::name =
    std::string("Qub<resize,double>");

template<>
const std::string qsl::Qubits<qsl::Type::Resize, float>::name =
    std::string("Qub<resize,float>");

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::Resize, Fp>::Qubits(unsigned nqubits_in) 
    : nqubits{ nqubits_in }, dim{ std::size_t(1) << nqubits },
      dim_max{ dim }, state(dim)
{
    // Make the all-zero state
    reset();   
}

template<std::floating_point Fp>
qsl::Qubits<qsl::Type::Resize, Fp>::Qubits(const std::vector<qsl::complex<Fp>> & state)
    : nqubits{ qsl::checkStateSize(state) }, 
      dim{ state.size() }, dim_max{ dim }, state{state}
{
    //std::cout << "nqubits = " << nqubits << std::endl;
    //std::cout << "dim = " << dim << std::endl;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::reset()
{
    for (std::size_t n=0; n<dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    state[0].real = 1;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::setState(const std::vector<qsl::complex<Fp>> & state_in)
{
    // For the resizeable simulator, you can set a different sized
    // state vector without causing an exception (unless the state is
    // an invalid size)
    dim = state_in.size();
    nqubits = checkStateSize(state_in);
   
    // Set the new state vector
    state = state_in;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::setBasisState(std::size_t index)
{
    // Clear the state - set all amplitudes to zero
    for (std::size_t n = 0; n < dim; n++) {
	state[n].real = 0;
	state[n].imag = 0;
    }
    // Set the amplitude for index to 1
    state[index].real = 1;
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::reallocateState()
{
    // This function should trim the state vector down to
    // the smallest possible size (i.e. set dim_max = dim)
    state.resize(dim);
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::appendQubit()
{
    addQubit(getNumQubits());
}

template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::addQubit(unsigned targ)
{
    // Check if there is any space left in the
    // state vector to allocate more qubits
    if (dim == dim_max) {
	// In this case, you need to make more space at the end
	// of the state vector. Call a resize operation which
	// address initialised memory
	state.resize(2*dim);
	// Update the max dim to the new state vector size
	dim_max = state.size();
    }

    // Now you need to shuffle the amplitudes so that they are
    // in the correct positions for the new qubit. This amounts
    // to doing the opposite of the operation in collapseOut.
    // You need to go backwards (otherwise you're going to
    // overwrite amplitudes with zeros... think about it...)
    std::size_t k = 1 << targ; // Stride length between blocks
    for (long long int n = (dim >> targ)-1; n >= 0 ; n--) {
	for (std::size_t p = 0; p < k; p++) {
	    // This is opposite to collapseOut. Use
	    // outcome = 0 to copy to the position of the
	    // zero amplitude.
	    std::size_t from = n*k + p;
	    std::size_t to = 2*n*k + p;
	    // This is the state to zero, corresponding to the
	    // one amplitude.
	    std::size_t to_zero = (2*n+1)*k + p;
	    state[to] = state[from];
	    state[to_zero].real = 0;
	    state[to_zero].imag = 0;
	}
    }    

    
    // Update the state vector dimension and number of qubits
    dim <<= 1;
    nqubits++;

}


/*
template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::operator = (const Qubits & old)
{
    // For the resizeable simulator, you can assign a different sized
    // simulator object without causing an exception (unless the state is
    // an invalid size)
    dim = old.dim;
    dim_max = old.dim_max;
    nqubits = old.nqubits;
    random = old.random;
    
    // Set the new state vector
    state = old.state;
}
*/
/// Get the state vector associated to the qubits
template<std::floating_point Fp>
std::vector<qsl::complex<Fp>> qsl::Qubits<qsl::Type::Resize, Fp>::getState() const
{
    auto first = std::begin(state);
    auto last = first + dim;
    return std::vector<qsl::complex<Fp>>(first, last);
}

template<std::floating_point Fp>
unsigned qsl::Qubits<qsl::Type::Resize, Fp>::getNumQubits() const
{
    return nqubits;
}


template<std::floating_point Fp>
void qsl::Qubits<qsl::Type::Resize, Fp>::print(std::ostream & os) const
{
    os << "Number of qubits = " << nqubits << std::endl;

    // Need to print manually because some of the state vector
    // at the end might be ignored
    for (std::size_t n = 0; n < dim; n++) {
	os << state[n] << std::endl;
    }
    os << std::endl;
}


// Explicit instantiations
template class qsl::Qubits<qsl::Type::Resize, float>;
template class qsl::Qubits<qsl::Type::Resize, double>;

template<> qsl::Random<float> qsl::Qubits<qsl::Type::Resize, float>::random{0,1};
template<> qsl::Random<double> qsl::Qubits<qsl::Type::Resize, double>::random{0,1};
