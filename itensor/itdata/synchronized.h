//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SYNCHRONIZED_H
#define __ITENSOR_SYNCHRONIZED_H

#include <mutex>

template<typename T>
class Synchronized
    {
    T data_;
    mutable std::mutex m_;
    public:

    using value_type = T;

    Synchronized() { }
    Synchronized(const value_type& data) : data_(data) { }
    Synchronized(value_type&& data) : data_(std::move(data)) { }

    Synchronized(const Synchronized& o) 
        { 
        set(o.data_);
        }

    Synchronized(Synchronized&& o) 
        { 
        set(std::move(o.data_));
        }

    Synchronized&
    operator=(const Synchronized& o) 
        { 
        set(o.data_);
        }

    Synchronized&
    operator=(Synchronized&& o) 
        { 
        set(std::move(o.data_));
        }

    explicit operator bool() const
        { 
        std::lock_guard<std::mutex> g(m_);
        return static_cast<bool>(data_);
        }

    value_type
    get() const 
        { 
        std::lock_guard<std::mutex> g(m_);
        return data_;
        }

    void
    set(const value_type& val)
        { 
        std::lock_guard<std::mutex> g(m_);
        data_ = val;
        }

    void
    set(value_type&& val)
        { 
        std::lock_guard<std::mutex> g(m_);
        data_ = std::move(val);
        }

    };

#endif
