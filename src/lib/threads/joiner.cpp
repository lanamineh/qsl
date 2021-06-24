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
 * \file joiner.cpp
 * \brief Implementation of the thread joiner
 *
 */

#include "qsl/threads/joiner.hpp"

JoinThreads::JoinThreads(std::vector<std::thread>& threads_in)
    : threads(threads_in) { }


JoinThreads::~JoinThreads()
{
    for(unsigned long i=0; i<threads.size(); ++i) {
	if(threads[i].joinable())
	    threads[i].join();
    }
}

