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
 * \file results.tpp
 * \brief Implementation of the Results class
 */

#include "qsl/benchmark/results.hpp"
#include <stdexcept>
#include <iostream>
#include <fstream>

namespace qsl {

    template<std::integral T, std::floating_point R>
    Results<T,R>::Results(const std::vector<std::string> & headings)
	: head{headings}
    { }

    template<std::integral T, std::floating_point R>
    template<typename S>
    std::vector<S> Results<T,R>::column(std::size_t index) const
    {
	std::vector<S> column;

	std::size_t rows = intcols.size();
    
	// Check for zero rows
	if (rows == 0) {
	    return column; // Return empty column
	}

	// If rows not zero, check whether int or fp column is requested
	std::size_t intnum = intcols[0].size();
	if (index < intnum) {
	    // Copy the integer column to columns
	    for (std::size_t r = 0; r < rows; r++) {
		column.push_back(intcols[r][index]);
	    }
	} else {
	    // Copy the floating point column to columns (correcting for index)
	    for (std::size_t r = 0; r < rows; r++) {
		column.push_back(fpcols[r][index - intnum]);
	    }
	}

	return column;
    }

    template<std::integral T, std::floating_point R>
    void Results<T, R>::addMeta(const std::string & line)
    {
	meta.push_back(line);
    }

    template<std::integral T, std::floating_point R>
    void Results<T, R>::addRow(const std::vector<T> & intcols_in,
			       const std::vector<R> & fpcols_in)
    {
	std::size_t nrows = intcols.size();

	if (nrows == 0) {
	    // Check total number of columns
	    if (intcols_in.size() + fpcols_in.size() != head.size()) {
		std::string msg = "Total number of columns does not match ";
		msg += "number of column headings in first call to addRow";
		throw std::logic_error(msg);
	    }
	} else {
	    // Check that the number of integer columns is correct
	    if (intcols[0].size() != intcols_in.size()) {
		std::string msg = "Incorrect number of integer columns in call ";
		msg += "to addRow";
		throw std::logic_error(msg);
	    }
	    // Check that the number of floating point columns is correct
	    if (fpcols[0].size() != fpcols_in.size()) {
		std::string msg = "Incorrect number of floating point columns in ";
		msg += "call to addRow";
		throw std::logic_error(msg);
	    }
	}
    
	// Store the columns
	intcols.push_back(intcols_in);
	fpcols.push_back(fpcols_in);
    }


    template<std::integral T, std::floating_point R>
    std::vector<std::string> Results<T,R>::headings() const
    {
	return head;
    }


    template<std::integral T, std::floating_point R>
    std::vector<std::string> Results<T, R>::getMeta() const
    {
	return meta;
    }

    template<std::integral T, std::floating_point R>
    template<typename S>
    S Results<T, R>::mean(const std::string & colname) const
    {
	// Look for colname index
	for (std::size_t index = 0; index < head.size(); index++) {
	    if (head[index] == colname) {
		return mean<S>(index);
	    }
	}

	// If you get here, the column name was not found
	throw std::logic_error("No such column name '" + colname + "' in Results"
			       " object");
    }

    template<std::integral T, std::floating_point R>
    template<typename S>
    S Results<T, R>::variance(const std::string & colname) const
    {
	// Look for colname index
	for (std::size_t index = 0; index < head.size(); index++) {
	    if (head[index] == colname) {
		return variance<S>(index);
	    }
	}

	// If you get here, the column name was not found
	throw std::logic_error("No such column name '" + colname + "' in Results"
			       " object");

    }

    template<std::integral T, std::floating_point R>
    template<typename S>
    S Results<T, R>::mean(std::size_t index) const
    {
	if (!(index < head.size())) {
	    std::string msg = "Index " + std::to_string(index) + " out of range ";
	    msg += "in call to mean";
	    throw std::out_of_range(msg);	
	}

	std::size_t nrows = intcols.size();
	S sum = 0;
	if (nrows != 0) {

	    // Get the number of integer and floating point columns
	    std::size_t nintcols = intcols[0].size();
	
	    if (index < nintcols) {
		for (std::size_t row = 0; row < nrows; row++) {
		    sum += intcols[row][index];
		}	    
	    } else {
		for (std::size_t row = 0; row < nrows; row++) {
		    sum += fpcols[row][index - nintcols];
		}
	    }
	}
	S average = sum/nrows;
	return average;
    }

    template<std::integral T, std::floating_point R>
    template<typename S>
    S Results<T,R>::variance(std::size_t index) const
    {
	if (!(index < head.size())) {
	    std::string msg = "Index " + std::to_string(index) + " out of range ";
	    msg += "in call to variance";
	    throw std::out_of_range(msg);	
	}
    
	std::size_t nrows = intcols.size();

	S x_bar = mean<S>(index); 

	S sum = 0;
	if (nrows != 0) {

	    // Get the number of integer and floating point columns
	    std::size_t nintcols = intcols[0].size();
	
	    if (index < nintcols) {
		for (std::size_t row = 0; row < nrows; row++) {
		    sum += std::pow(intcols[row][index] - x_bar, 2);
		}	    
	    } else {
		for (std::size_t row = 0; row < nrows; row++) {
		    sum += std::pow(fpcols[row][index - nintcols] - x_bar, 2);
		}
	    }
	}
	S variance = sum/(nrows - 1);
	return variance;
    }

    template<std::integral T, std::floating_point R>
    Results<T,R> Results<T,R>::subset(const std::vector<std::size_t> & indices) const
    {
	// Make the new list of headings
	std::vector<std::string> new_head;
	for (std::size_t n = 0; n < indices.size(); n++) {
	    if (indices[n] < head.size()) {
		std::size_t col = indices[n];
		new_head.push_back(head[col]);
	    } else {
		throw std::out_of_range("Index " + std::to_string(indices[n])
					+ " is out of range in call to subset");
	    }
	}
    
	// Make a new results object
	Results<T,R> results{new_head};

	// Make new intcols
	std::size_t nrows = intcols.size();

	if (nrows != 0) {
	    ///\todo Copy columns to new Results object
	    for (std::size_t row = 0; row < nrows; row++) {

		std::size_t intlim = intcols[0].size();
		std::size_t fplim = fpcols[0].size();

		// Extract row elements to form new sub rows for
		std::vector<T> intcols_row;
		std::vector<R> fpcols_row;
		for (std::size_t n = 0; n < indices.size(); n++) {

		    // Check if index points to integer or floating point column
		    if (indices[n] < intlim) {
			std::size_t col = indices[n];
			intcols_row.push_back(intcols[row][col]);
		    } else if ((indices[n] - intlim) < fplim) {
			std::size_t col = indices[n] - intlim;
			fpcols_row.push_back(fpcols[row][col]);		    
		    }		
		}

		// Push the new rows to the new results object 
		results.intcols.push_back(intcols_row);
		results.fpcols.push_back(fpcols_row);
	    }
	}
	    
	return results;
    
    }


    template<std::integral T, std::floating_point R>
    void Results<T, R>::print() const
    {
	std::cout << std::endl;

	// Print the metadata
	std::cout << "--- META-DATA ---" << std::endl;
	for(std::size_t n = 0; n < meta.size(); n++) {
	    std::cout << meta[n] << std::endl;
	}
    
	// Print the columns
	std::cout << std::endl << "--- TABLE BODY ---" << std::endl;
	for(std::size_t n = 0; n < head.size()-1; n++) {
	    std::cout << head[n] << ", ";
	}
	// Final column headings (no trailing comma)
	std::cout << head.back() << std::endl;

	// Print the table body
	std::size_t nrows = intcols.size();
	if (nrows != 0 ) {
	    // Get the number of integer and floating point columns
	    std::size_t nintcols = intcols[0].size();
	    std::size_t nfpcols = fpcols[0].size();
	
	    // Print all rows
	    for(std::size_t row = 0; row < nrows; row++) {
		// Loop over integer columns
		for (std::size_t col = 0; col < nintcols; col++) {
		    std::cout << intcols[row][col] << ", ";
		}
		// Loop over floating point columns
		for (std::size_t col = 0; col < nfpcols - 1; col++) {
		    std::cout << fpcols[row][col] << ", ";
		}
		// Print last column
		std::cout << fpcols[row].back() << std::endl;
	    }
	}
	std::cout << "------------------" << std::endl;

    }

    template<std::integral T, std::floating_point R>
    void Results<T, R>::writeToFile(const std::string & filename) const
    {
	std::ofstream file;
	file.open(filename);
    
	// Print the metadata
	for(std::size_t n = 0; n < meta.size(); n++) {
	    file << "# " << meta[n] << std::endl;
	}
    
	// Print the columns
	file << std::endl;    
	for(std::size_t n = 0; n < head.size()-1; n++) {
	    file << head[n] << " ";
	}
	// Final column headings (no trailing comma)
	file << head.back() << std::endl;

	// Print the table body
	std::size_t nrows = intcols.size();
	if (nrows != 0 ) {
	    // Get the number of integer and floating point columns
	    std::size_t nintcols = intcols[0].size();
	    std::size_t nfpcols = fpcols[0].size();
	
	    // Print all rows
	    for(std::size_t row = 0; row < nrows; row++) {
		// Loop over integer columns
		for (std::size_t col = 0; col < nintcols; col++) {
		    file << intcols[row][col] << " ";
		}
		// Loop over floating point columns
		for (std::size_t col = 0; col < nfpcols - 1; col++) {
		    file << fpcols[row][col] << " ";
		}
		// Print last column
		file << fpcols[row].back() << std::endl;
	    }
	}
    }
    
}
