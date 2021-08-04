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
 * \file timer.cpp
 * \brief Contains functions from timer.hpp
 */

#include "qsl/utils/timer.hpp"
#include <ctime>

void qsl::Timer::start()
{
    time1 = clock_t::now();
}

void qsl::Timer::stop()
{
    time2 = clock_t::now();
    elapsed = std::chrono::duration_cast<second_t>(time2 - time1).count();
}

double qsl::Timer::getElapsed()
{
    return elapsed;
}

std::string qsl::Timer::printElapsed()
{
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1);
    auto h = std::chrono::duration_cast<std::chrono::hours>(ms);
    ms -= h;
    auto m = std::chrono::duration_cast<std::chrono::minutes>(ms);
    ms -= m;
    auto s = std::chrono::duration_cast<std::chrono::seconds>(ms);
    ms -= s;
    
    std::string timeTaken;
    timeTaken = std::to_string(h.count()) + "h "
	+ std::to_string(m.count()) + "m "
	+ std::to_string(s.count()) + "s "
	+ std::to_string(ms.count()) + "ms";

    return timeTaken;
}

std::string qsl::Timer::getCurrentTime()
{
    auto timenow = clock_t::to_time_t(clock_t::now());
    return std::string(ctime(&timenow));
}

