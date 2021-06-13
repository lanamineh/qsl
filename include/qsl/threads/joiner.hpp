/**
 * \file joiner.hpp
 * \brief A thread joiner
 *
 * This class is taken from 'C++ Concurrency in Action'
 */

#ifndef JOINER_HPP
#define JOINER_HPP

#include "queue.hpp"
#include <thread>

class JoinThreads
{
    std::vector<std::thread>& threads;
public:
    explicit JoinThreads(std::vector<std::thread>& threads_in);
    ~JoinThreads();
};

#endif
