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
 * \file quantum.cpp
 * \brief Implementation of functions in quantum.hpp
 *
 */

#define _USE_MATH_DEFINES // For MSVC, to use M_PI

#include <random>
#include <vector>
#include <string>
#include "qsl/utils/complex.hpp"
//#include <cmath>

template<typename Fp = double>
complex<Fp> innerProduct(const std::vector<complex<Fp>> & v,
			 const std::vector<complex<Fp>> & w)
{
    if(v.size() != w.size()) {
	std::string msg = "Cannot compute inner product between vectors ";
	msg += "of different sizes: ";
	msg += std::to_string(v.size()) + " and " + std::to_string(w.size());
	throw std::logic_error(msg);
    }

    complex<Fp> result{ 0,0 };
    for(std::size_t n=0; n<v.size(); n++) {
	result.real += v[n].real * w[n].real + v[n].imag * w[n].imag;
	result.imag += v[n].real * w[n].imag - v[n].imag * w[n].real;
    }

    return result;
}

template<typename Fp = double>
Fp norm(const std::vector<complex<Fp>> & v)
{
    // Find the norm of the vector
    Fp inner_prod = 0;
    for (std::size_t i = 0; i < v.size(); i++) {
	inner_prod += v[i].real * v[i].real + v[i].imag * v[i].imag;
    }

    return std::sqrt(inner_prod);
    
}

/**
 * \brief Compute arccos after checking argument
 *
 * The function is necessary in case a number rounds 
 * to something bigger than 1 or smaller than -1
 *
 */
double acosSafe(double value) {
    if (value <= -1.0) {
        return M_PI;
    } else if (value >= 1.0) {
        return 0;
    } else {
        return std::acos(value);
    }
}

template<typename Fp = double>
Fp fubiniStudy(const std::vector<complex<Fp>> & v,
	       const std::vector<complex<Fp>> & w)
{
    Fp numerator = abs(innerProduct(v,w));
    Fp denominator = norm(v) * norm(w); 
    return acosSafe(numerator/denominator);
}

template<typename Fp = double>
Fp normalise(std::vector<complex<Fp>> &state)
{
    // Find the norm of the vector
    Fp factor = norm(state);
   
    // Divide by the norm;
    for (std::size_t i = 0; i < state.size(); i++) {
	state[i].real /= factor;
	state[i].imag /= factor;
    }

    return factor;
}

/// Explicit instantiations
template float fubiniStudy(const std::vector<complex<float>> & v,
			   const std::vector<complex<float>> & w);
template double fubiniStudy(const std::vector<complex<double>> & v,
			    const std::vector<complex<double>> & w);

template float normalise(std::vector<complex<float>> &state);
template double normalise(std::vector<complex<double>> &state);
