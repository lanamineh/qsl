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

