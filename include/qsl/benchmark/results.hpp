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
 * \file results.hpp
 * \brief Contains a simple class for storing a table of results with metadata
 */

#ifndef RESULTS_HPP
#define RESULTS_HPP

#include <string>
#include <map>
#include <vector>

namespace qsl {
    /**
     * \brief A simple results class for holding a table of columns and metadata 
     *
     * A simple class which can a table or numerical data. Each row is a 
     * std::tuple, meaning that the columns can be different numerical types 
     * (e.g. double or int). It is necessary that the columns types be specified 
     * up front as template arguments to the Results type when the object is 
     * constructed. The corresponding column heading names must be specified in
     * the constructor. After construction, rows can be added to the table by 
     * passing a std::tuple to the addRow member function.
     *
     * The class also stores metadata strings, using the addMeta member function.
     * The table along with the metadata can be written to a file using 
     * writeToFile.
     *
     */
    template<std::integral T, std::floating_point R>
    class Results
    {
	std::vector<std::string> meta; ///< Name of gate, sims being compared, etc.
	std::vector<std::string> head; ///< Column headings
	std::vector<std::vector<T>> intcols; ///< Integer data columns
	std::vector<std::vector<R>> fpcols; ///< Floating point data columns
    
    public:

	/**
	 * \brief Make an object with the given column headings
	 * \param headings The columns headings, in order, begining with the first
	 */
	Results(const std::vector<std::string> & headings);	
    
	/**
	 * \brief Get a single column as a std::vector
	 *
	 * This function return a single column of data from the table. 
	 * The template parameter sets the type of the returned standard
	 * vector. It is important to choose a type which is compatible
	 * with the column type, otherwise compiler warnings or errors
	 * may result.
	 *
	 * \param index The column index to return
	 */
	template<typename S>
	std::vector<S> column(std::size_t index) const;

	/**
	 * \brief Get a single column as a std::vector
	 *
	 * This function gets a single column of data from the object,
	 * based on the column name. A std::logic_error will be thrown
	 * if the column name is not known.
	 *
	 * \param index The column index to return
	 */
	template<typename S>
	std::vector<S> column(const std::string & colname) const;
    
	/**
	 * \brief Get a subset of the columns as a new Results object
	 *
	 * This function returns a subset of the columns in a new 
	 * Results object with the same integer and floating point types.
	 * The list of indices must be within range, otherwise a 
	 * std::out_of_range exception is thrown.
	 *
	 * The return object preserves the column headings of the subset,
	 * but no metadata is included in the returned object.
	 *
	 * \param indices The list of column indices to include in the returned
	 * Results object.
	 * \return A new Results object containing the subset of columns.
	 */
	Results<T,R> subset(const std::vector<std::size_t> & indices) const;

	/**
	 * \brief Get a subset of the columns as a new Results object
	 *
	 * This function is the same as the other subset member function,
	 * but gets columns by name instead of by index. A std::logic_error
	 * is thrown if any of the column names are not recognised.
	 *
	 * \param colnames The list of column names to include in the returned
	 * Results object.
	 * \return A new Results object containing the subset of columns.
	 */
	Results<T,R> subset(const std::vector<std::string> & colnames) const;
    
	/**
	 * \brief Get the column headings as a standard vector of strings
	 * \return A standard vector of strings containing the column names.
	 */
	std::vector<std::string> headings() const;
    
	/**
	 * \brief Compute a column mean
	 *
	 * Compute and return the average value of a column of data. The
	 * column is specified using the index parameter. The return
	 * type of the data defaults is set using the template parameter. This
	 * might be useful for average a column of ints and returning the 
	 * result as a double. It is important to make sure that the return
	 * type is compatible with the column type, otherwise compiler warnings
	 * or errors may result.
	 *
	 * \param index The column index to average
	 * \param S The type for the returned average
	 * \return The average value of the selected column
	 */
	template<typename S>
	S mean(std::size_t index) const;

	/**
	 * \brief Compute a column mean
	 *
	 * Compute and return the average value of a column of data. The
	 * column is specified using the name parameter. See the other
	 * mean member function for information about the template parameter.
	 * The function will throw std::out_of_range if the column name is
	 * not recognised.
	 *
	 * \param colname The name of the column to average
	 * \param S The type for the returned average
	 * \return The average value of the selected column
	 */    
	template<typename S>
	S mean(const std::string & colname) const;

	/**
	 * \brief Compute an unbiased esimate of a column variance
	 *
	 * Compute and return the variance of a column of data. The
	 * column is specified using the index parameter. The return
	 * type behaves the same way as for the mean(std::size_t) member 
	 * function.
	 *
	 * The function computes the unbiased variance using the formula
	 *
	 * \f[
	 * \text{Var} = \frac{1}{n-1} \sum_{i=0}^n (x_i - \overline{x})^2
	 * \f]
	 * 
	 * where \f$ {x_i} \f$ are the column data, \f$ n \f$ is the column 
	 * length and \f$\overline{x}\f$ is the column mean
	 *
	 * \param index The column index for the variance
	 * \param S The type for the returned variance.
	 * \return The variance of the select column
	 *
	 */
	template<typename S>
	S variance(std::size_t index) const;

	/**
	 * \brief Compute an unbiased esimate of a column variance
	 *
	 * Compute and return the variance of a column of data, 
	 * specifying the column name using the colname parameter.
	 * See the other variance(std::size_t) member function for
	 * more information about the template parameter and the
	 * calculation of the variance. The function will throw
	 * std::out_of_range if the column name is not recognised.
	 *
	 * \param colname The column name for the variance
	 * \param S The type for the returned variance.
	 * \return The variance of the select column
	 *
	 */
	template<typename S>
	S variance(const std::string & colname) const;    
    
	/**
	 * \brief Add a whole row of data
	 *
	 * Add a full row of data to the Results object, specified as
	 * a standard vector of integer data and a standard vector
	 * of floating point data. The types of the columns should match
	 * that of the class. 
	 *
	 * The first call to the function sets the number of integer
	 * and floating point columns. The total number of columns must
	 * add up to the number of headings, otherwise a std::logic_error
	 * will be thrown. For subsequent calls, the same lengths for
	 * integer and floating point columns must be used, otherwise
	 * std::logic_error will be thrown.
	 *
	 * \param intcols The integer data part of the row
	 * \param fpcols The floating point data part of the row
	 */
	void addRow(const std::vector<T> & intcols, const std::vector<R> & fpcols);
    
	/**
	 * \brief Add a line of meta data 
	 */
	void addMeta(const std::string & line);

	/**
	 * \brief Get the metadata lines 
	 */
	std::vector<std::string> getMeta() const;
    
	/**
	 * \brief Write the results to a file
	 */
	void writeToFile(const std::string & filename) const;

	/**
	 * \brief Print the results to the screen
	 */
	void print() const;
    };
}
#include "qsl/benchmark/results.tpp"

#endif 
