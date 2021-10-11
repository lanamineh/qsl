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
 * \file pool.cpp
 * \brief Implementation of the thread pool
 *
 */

#include "qsl/threads/pool.hpp"

void Pool::worker_thread()
{
    while(!done) {
	FunctionWrapper task;
	if(work_queue.try_pop(task)) {
	    task();
	}
	else {
	    std::this_thread::yield();
	}
    }
}

// void Pool::worker_thread()
// {
//     while(!done) {
// 	std::function<void()> task;
// 	if(work_queue.try_pop(task)) {
// 	    task();
// 	}
// 	else {
// 	    std::this_thread::yield();
// 	}
//     }
// }

Pool::Pool()
    : done(false), joiner(threads)
{
    unsigned const thread_count=std::thread::hardware_concurrency();
    std::cout << "Creating thread pool with " << thread_count
	      << " threads" << std::endl;
    try {
	for(unsigned i = 0; i < thread_count; ++i) {
	    threads.push_back(
		std::thread(&Pool::worker_thread,this));
	}
    } catch(...) {
	done=true;
	    throw;
    }
}

Pool::~Pool()
{
    done=true;
}
