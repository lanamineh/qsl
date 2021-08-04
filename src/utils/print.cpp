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
