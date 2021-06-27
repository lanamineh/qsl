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
 * \file complex.hpp
 * \brief A complex type that bahaves a bit like std::complex
 *
 * 
 */

#ifndef COMPLEX_HPP
#define COMPLEX_HPP

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
