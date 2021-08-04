/*
 *  Copyright 2021 Lana Mineh and John Scott
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
