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
 * \file misc.cpp
 * \brief Implementation for functions in misc.hpp
 *
 */

#include "qsl/utils/misc.hpp"
#include <bitset>
#include <string>
#include <iostream>
#include <iomanip>
#include "qsl/utils/colours.hpp"

namespace qsl {

    // Use __builtin_ctz with gcc/clang/icc
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define ctz __builtin_ctz

    // Catch other cases
#else
#error The compiler you are using is not recognised
#endif


    void next(std::size_t & x) {

	// t gets x's least significant 0 bits set to 1
	std::size_t t = x | (x - 1);

	// Next set to 1 the most significant bit to change,
	// set to 0 the least significant ones, and add the necessary 1 bits.
	x = (t + 1) | (((~t & -~t) - 1) >> (ctz(x) + 1));  
    }

    unsigned hammingWeight(std::size_t n)
    {
	unsigned count = 0; 
	while (n) { 
	    count += n & 1; 
	    n >>= 1; 
	} 
	return count; 
    }


    unsigned choose(unsigned n, unsigned k)
    {
	unsigned ans = 1;
	for (unsigned i = 1; i <= k; i++) {
	    ans *= (n - k + i);
	    ans /= i;
	}
	return ans;
    }

}
