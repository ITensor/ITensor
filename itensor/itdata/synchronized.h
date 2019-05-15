//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
