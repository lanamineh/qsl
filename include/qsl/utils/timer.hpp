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
 * \file timer.hpp
 * \brief Header file for timing functions.
 */

#ifndef QSL_TIMER_HPP
#define QSL_TIMER_HPP

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

