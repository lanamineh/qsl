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
 * \file print.cpp
 * \brief Implementation of printing related functions
 *
 */

#include "qsl/utils/misc.hpp"
#include <bitset>
#include <string>
#include <iostream>
#include <iomanip>
#include "qsl/utils/colours.hpp"

namespace qsl {

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
	std::cout << "0b" << std::bitset<10>(x) << std::endl;
    }

    
}


template<std::floating_point Fp>
std::ostream & operator<< (std::ostream & stream, qsl::complex<Fp> val)
{
    stream << "(" << val.real << ", " << val.imag << ")";
    return stream;
}


template std::ostream & operator<< (std::ostream & stream,
				    qsl::complex<float> val);
template std::ostream & operator<< (std::ostream & stream,
				    qsl::complex<double> val);
