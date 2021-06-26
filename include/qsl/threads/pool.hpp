/**
 * \file pool.hpp
 * \brief Contains a thread pool for use by the simulator
 *
 * This class is taken from 'C++ Concurrency in Action'
 */

#ifndef POOL_HPP
#define POOL_HPP

#include "queue.hpp"
#include "joiner.hpp"

#include <vector>
#include <atomic>
#include <thread>
#include <functional>
#include <iostream>
#include <future>

class FunctionWrapper
{
    struct impl_base {
	virtual void call()=0;
	virtual ~impl_base() {}
    };
    std::unique_ptr<impl_base> impl;
    template<typename F>
    struct impl_type : impl_base
    {
	F f;
	impl_type(F&& f_): f(std::move(f_)) {}
	void call() { f(); }
    };
public:
    template<typename F>
    FunctionWrapper(F&& f):
	impl(new impl_type<F>(std::move(f)))
	{}
    void operator()() { impl->call(); }
    FunctionWrapper() = default;
    FunctionWrapper(FunctionWrapper&& other):
	impl(std::move(other.impl))
	{}

    FunctionWrapper& operator=(FunctionWrapper&& other)
	{
	    impl=std::move(other.impl);
	    return *this;
	}
    FunctionWrapper(const FunctionWrapper&)=delete;
    FunctionWrapper(FunctionWrapper&)=delete;
    FunctionWrapper& operator=(const FunctionWrapper&)=delete;
};

class Pool
{
    std::atomic_bool done;
    ThreadsafeQueue<FunctionWrapper> work_queue;
    std::vector<std::thread> threads;
    JoinThreads joiner;
    void worker_thread();
public:
    Pool();
    ~Pool();
    
    template<typename FunctionType>
    std::future<typename std::result_of<FunctionType()>::type>
    submit(FunctionType f) {
	
	using result_type = typename std::result_of<FunctionType()>::type;
	std::packaged_task<result_type()> task(std::move(f));
	std::future<result_type> res(task.get_future());
	work_queue.push(std::move(task));
	return res;
    }
    
};

#endif
