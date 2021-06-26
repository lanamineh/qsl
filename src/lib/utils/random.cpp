/* 
 * Copyright (C) 2020 Lana Mineh and John Scott.
 *
 * This file is part of QSL, the quantum computer simulator.
 *
 * QSL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QSL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QSL.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * \file random.cpp
 * \brief Implementation of random number related functions
 *
 *
 */

#include "qsl/utils/random.hpp"
#include "qsl/utils/quantum.hpp"
#include "qsl/utils/misc.hpp"

namespace qsl {
 
    template<typename Fp>
    Random<Fp>::Random(Fp a, Fp b)
    {
	///\todo Is std::random_device OK?
	std::random_device r;
	
	generator = std::mt19937(r());
	distribution = std::uniform_real_distribution<Fp>(a,b);
    }

    template<typename Fp>
    Fp Random<Fp>::getNum() {
	return distribution(generator);
    }


template<typename Fp>
std::vector<complex<Fp>> makeRandomState(std::uint8_t nqubits)
{
    std::size_t dim = 1 << nqubits;

    qsl::Random<Fp> random(-1,1);

    
    
    std::vector<complex<Fp>> state;
    for(std::size_t i=0; i<dim; i++) {
	Fp val_real = random.getNum();
	Fp val_imag = random.getNum();
	state.push_back(complex(val_real, val_imag));
    }

    // Normalise the state vector
    normalise(state);
    
    return state;
}

///\todo This needs testing
template<typename Fp>
std::vector<complex<Fp>> makeRandomNPState(std::uint8_t nqubits)
{
    std::size_t dim = 1 << nqubits;

    qsl::Random<Fp> random(-1,1);

    // Make a random number of ones
    std::random_device r;    
    std::mt19937 generator{r()};
    std::uniform_int_distribution<unsigned> distribution(1,nqubits-1);
    unsigned nones = nqubits/2;//distribution(generator);

    std::size_t x = (1ULL << nones) - 1;
    std::size_t end = x << (nqubits - nones);

    // Make state vector with correct length
    std::vector<complex<Fp>> state(dim);
    
    while (x <= end) {

	Fp val_real = random.getNum();
	Fp val_imag = random.getNum();
	state[x] = complex(val_real, val_imag);	
	next(x);
    }
    
    // Normalise the state vector
    normalise(state);
    
    return state;
}


template<typename Fp>
std::vector<Fp> makeRandomPhases(int test_len)
{
    std::vector<Fp> phase_list;

    qsl::Random<Fp> random(-M_PI, M_PI);
    
    for (int i = 0; i < test_len; i++) {
	phase_list.push_back(random.getNum());
    }
    return phase_list;
}

/// Explicit instantiations
template class qsl::Random<float>;
template class qsl::Random<double>;

template std::vector<complex<float>> makeRandomState(std::uint8_t nqubits);
template std::vector<complex<double>> makeRandomState(std::uint8_t nqubits);

template std::vector<complex<float>> makeRandomNPState(std::uint8_t nqubits);
template std::vector<complex<double>> makeRandomNPState(std::uint8_t nqubits);

template std::vector<float> makeRandomPhases(int nqubits);
template std::vector<double> makeRandomPhases(int nqubits);

}
