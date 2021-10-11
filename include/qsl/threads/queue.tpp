/**
 * \file queue.tpp
 * \brief Implementation of thread safe queue
 *
 */

template<typename T>
ThreadsafeQueue<T>::ThreadsafeQueue() {}

template<typename T>
ThreadsafeQueue<T>::ThreadsafeQueue(ThreadsafeQueue const& other)
{
    std::lock_guard<std::mutex> lk(other.mut);
    data_queue=other.data_queue;
}

template<typename T>
void ThreadsafeQueue<T>::push(T && new_value)
{
    std::lock_guard<std::mutex> lk(mut);
    data_queue.push(std::move(new_value));
}

template<typename T>
void ThreadsafeQueue<T>::wait_and_pop(T& value)
{
    std::unique_lock<std::mutex> lk(mut);
    data_cond.wait(lk,[this]{return !data_queue.empty();});
    value=data_queue.front();
    data_queue.pop();
}

template<typename T>
std::shared_ptr<T> ThreadsafeQueue<T>::wait_and_pop()
{
    std::unique_lock<std::mutex> lk(mut);
    data_cond.wait(lk,[this]{return !data_queue.empty();});
    std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
    data_queue.pop();
    return res;
}

template<typename T>
bool ThreadsafeQueue<T>::try_pop(T& value)
{
    std::lock_guard<std::mutex> lk(mut);
    if(data_queue.empty())
	return false;
    value=std::move(data_queue.front());
    data_queue.pop();
    return true;
}

template<typename T>
std::shared_ptr<T> ThreadsafeQueue<T>::try_pop()
{
    std::lock_guard<std::mutex> lk(mut);
    if(data_queue.empty())
	return std::shared_ptr<T>();
    std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
    data_queue.pop();
    return res;
}

template<typename T>
bool ThreadsafeQueue<T>::empty() const
{
    std::lock_guard<std::mutex> lk(mut);
    return data_queue.empty();
}
