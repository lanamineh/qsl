/**
 * \file queue.hpp
 * \brief Contains a threadsafe queue
 *
 * This class is taken from 'C++ Concurrency in Action'
 */

#ifndef QUEUE_HPP
#define QUEUE_HPP

#include <mutex>
#include <queue>
#include <condition_variable>
#include <memory>

template<typename T>
class ThreadsafeQueue
{
private:
    mutable std::mutex mut;
    std::queue<T> data_queue;
    std::condition_variable data_cond;
public:
    ThreadsafeQueue();
    ThreadsafeQueue(ThreadsafeQueue const& other);
    void push(T && new_value);
    void wait_and_pop(T& value);
    std::shared_ptr<T> wait_and_pop();
    bool try_pop(T& value);
    std::shared_ptr<T> try_pop();
    bool empty() const;
};

#include "queue.tpp"

#endif
