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

// Make a windows version of __builtin_ctz for the
// next function
#ifdef _MSC_VER
#include <intrin.h>
#include <windows.h>
uint32_t __inline ctz(uint32_t value)
{
    DWORD trailing_zero = 0;
    if (_BitScanForward(&trailing_zero, value)) {
        return trailing_zero;
    } else {
        // This is undefined for __builtin_ctz
        return 32;
    }
}

// Use __builtin_ctz with gcc/clang/icc
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
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

void printRow(std::string col1, double time1, double time2)
{   
    double diff = time2 - time1;

    // Compute the percentage slower or faster  
    double percent = 100*time1/time2;

    // Save iostream format flags
    std::ios_base::fmtflags f( std::cout.flags() );
    
    // Print results
    std::cout << col1 << ", "
	      << std::fixed << std::setprecision(2)
	      << 1000*time1 << "ms, " << 1000*time2 << "ms, ";
    if(diff > 0) {
	std::cout << Colour::GREEN;
    } else {
	std::cout << Colour::YELLOW;
    }
    std::cout << 1000*diff << std::fixed << "ms (" << percent << "%)"
       << Colour::RESET << std::endl;

    // Restore iostream flags
    std::cout.flags(f);
}

void printRow(double time1, double time2)
{   
    double diff = time2 - time1;

    // Compute the percentage slower or faster  
    double percent = 100*time1/time2;

    // Save iostream format flags
    std::ios_base::fmtflags f( std::cout.flags() );
    
    // Print results
    std::cout << std::fixed << std::setprecision(2)
	      << 1000*time1 << "ms, " << 1000*time2 << "ms, ";

    if(diff > 0) {
	std::cout << Colour::GREEN;
    } else {
	std::cout << Colour::YELLOW;
    }
    std::cout << 1000*diff << std::fixed << "ms (" << percent << "%)"
	      << Colour::RESET << std::endl;

    // Restore iostream flags
    std::cout.flags(f);
}


void printBin(const std::string & name, unsigned int x) {
    std::cout << name << " = ";
    std::cout << "0b" << std::bitset<5>(x) << std::endl;
}

template<typename Fp>
unsigned checkStateSize(const std::vector<complex<Fp>> & state) {

    // Counts the number of ones in the binary representation of dim
    // If this number is not 1, then dim is not a power of 2 and
    // is therefore invalid
    unsigned one_count = 0;
    // Counts the length of the binary representation of dim
    unsigned len = 0;
    
    std::size_t n = state.size();
    while(n > 0) {
	len += 1;
	one_count += n & 1;
	n >>= 1;
    }
    
    if(one_count != 1) {
	throw std::logic_error("State vector size is not a power of two");
    }

    return len - 1;

}


// Explicit instantiations
template unsigned checkStateSize(const std::vector<complex<float>> & state);
template unsigned checkStateSize(const std::vector<complex<double>> & state);


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
