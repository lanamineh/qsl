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

    std::vector<std::size_t> generateLookup(unsigned len, unsigned ones)
    {
	std::vector<std::size_t> list;
	// Starting number (ones in the least significant positions)
	std::size_t x = (1ULL << ones) - 1ULL;
	// Last number (ones in the most significant positions)
	std::size_t end = x << (len - ones);
	
	// Loop through all the other numbers
	while(x <= end) {
	    list.push_back(x);
	    qsl::next(x);
	}
	
	return list;
    }


    
}
