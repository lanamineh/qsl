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
 * \file timer.hpp
 * \brief Header file for timing functions.
 */

#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <string>

namespace qsl {

/**
 * \brief A class for timing things.
 */
    class Timer
    {
    private:
	// Type aliases to make accessing nested type easier
	/// Typedef for system clock
	using clock_t = std::chrono::system_clock;
	/// Typedef for duration 
	using second_t = std::chrono::duration<double, std::ratio<1>>; 
  
	std::chrono::time_point<clock_t> time1;  ///< Start time
	std::chrono::time_point<clock_t> time2;  ///< End time
	double elapsed;  ///< Elapsed time in seconds

    public:
	/**
	 * \brief Class constructor, sets elapsed time to 0.
	 */
	Timer() : elapsed(0) {}

	/**
	 * \brief Start the timer.
	 */
	void start();

	/**
	 * \brief Stop the timer.
	 */
	void stop();
    
	/**
	 * \brief Access function, get the elapsed time in seconds. 
	 * Note that the timer must be stopped.
	 */
	double getElapsed();
    
	/**
	 * \brief Function to return the elapsed time in terms of hours, 
	 * minutes, seconds and milliseconds.
	 */
	std::string printElapsed();
    
	/**
	 * \brief Function to return current date and time.
	 */
	static std::string getCurrentTime();
    };

}
#endif

