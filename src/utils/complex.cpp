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
 * \file complex.cpp
 * \brief Implementation of functions relating to the complex type
 *
 */

#include "qsl/utils/complex.hpp"
#include <cmath>

namespace qsl {

    template<std::floating_point Fp>
    complex<Fp>::complex(Fp r, Fp i)
    {
	real = r;
	imag = i;
    }

    template<std::floating_point Fp>
    complex<Fp>::complex() : complex(0,0) { }

    template<std::floating_point Fp>
    complex<Fp> operator - (const complex<Fp> & a, const complex<Fp> & b)
    {
	return complex{ a.real - b.real, a.imag - b.imag};
    }

    template<std::floating_point Fp>
    Fp abs(const complex<Fp> & a)
    {
	return std::sqrt(a.real*a.real + a.imag*a.imag);
    }

    /// Explicit template instantiations
    template struct complex<float>;
    template struct complex<double>;

    template float abs(const complex<float> & a);
    template double abs(const complex<double> & a);

    template complex<float> operator - (const complex<float> & a,
					const complex<float> & b);
    template complex<double> operator - (const complex<double> & a,
					 const complex<double> & b);
	
    
}
