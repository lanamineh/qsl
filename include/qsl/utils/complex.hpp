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
 * \file complex.hpp
 * \brief A complex type that bahaves a bit like std::complex
 *
 * 
 */

#ifndef QSL_COMPLEX_HPP
#define QSL_COMPLEX_HPP

#include <ostream>

namespace qsl {

    /**
     * \brief Struct containing a single complex number
     *
     * This class is deliberately minimal to avoid performance problems.
     *
     */
    template<std::floating_point Fp = double>
    struct complex
    {
	Fp real;
	Fp imag;

	/// Construct a comple number from its real and imaginary parts
	complex(Fp r, Fp i);
	complex(); ///< Make the complex number zero
    };

    // Deduction guide for complex constructor
    template<std::floating_point Fp>
    complex(Fp, Fp) -> complex<Fp>;

    /**
     * \brief Compute the absolute value of a complex number
     *
     * It is important this has the same name as std::abs so
     * that it can be used in template functions that work
     * for both double and complex types.
     *
     * \todo Maybe it would be good to add an absSquared?
     * 
     */
    template<std::floating_point Fp = double>
    Fp abs(const complex<Fp> & a);

    /**
     * \brief Compute the difference between two complex numbers
     *
     * This function is not meant to be fast
     *
     */
    template<std::floating_point Fp = double>
    complex<Fp> operator - (const complex<Fp> & a, const complex<Fp> & b);

}

/**
 * \brief Print a complex number to an output stream
 */
template<std::floating_point Fp = double>
std::ostream & operator<< (std::ostream & stream, qsl::complex<Fp> val);


#endif
