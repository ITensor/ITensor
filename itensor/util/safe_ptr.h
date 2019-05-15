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
#ifndef __ITENSOR_SAFE_PTR_H
#define __ITENSOR_SAFE_PTR_H

#include "itensor/util/print.h"
#include "itensor/util/stdx.h"

namespace itensor {

template<typename T>
class SafePtr
    {
    public:
    using value_type = stdx::remove_const_t<T>;
    using pointer = T*;
    using reference = T&;
    using const_pointer = const T*;
    private:
    pointer p_ = nullptr;
    size_t offset_ = 0;
    size_t offset_end_ = 0;
    public:

    //Defined non-explicit so we can write
    //SAFE_PTR_OF(T) sp = nullptr;
    SafePtr(T* pt = nullptr)
        {
        if(pt) throw std::runtime_error("SafePtr: single-arg constructor only accepts nullptr");
        }

    SafePtr(T* pt, size_t offset_end)
        : p_(pt), offset_(0), offset_end_(offset_end)
        { 
        //if(!p_) throw std::runtime_error("SafePtr: pointer is null");
        }

    SafePtr(T* pt, size_t offset, size_t offset_end)
        : p_(pt), offset_(offset), offset_end_(offset_end)
        { 
        if(!p_) throw std::runtime_error("SafePtr: pointer is null");
        }

    //Allow automatic conversion SafePtr<T> -> SafePtr<const T>
    SafePtr(SafePtr<value_type> const& P)
        { 
        operator=(P);
        }
    SafePtr&
    operator=(SafePtr<value_type> const& P)
        {
        p_ = P.getStart();
        offset_ = P.offset();
        offset_end_ = P.offsetEnd();
        return *this;
        }

    size_t
    offset() const { return offset_; }
    size_t
    offsetEnd() const { return offset_end_; }
    bool
    validOffset() const { return (offset_ < offset_end_); }
    size_t
    range() const 
        { 
        if(validOffset()) return (offset_end_-offset_);
        return 0ul;
        }

    pointer
    get() const 
        { 
        if(!p_) return nullptr;
        return (p_+offset_); 
        }

    pointer
    safeGet(size_t expected_range)
        {
        if(!p_) 
            {
            throw std::runtime_error("SafePtr: dereferencing null pointer");
            }
        if(!validOffset())
            {
            auto error_msg = 
            format("SafePtr: offset >= offset_end (%d >= %d)",
                   offset_,offset_end_);
            throw std::runtime_error(error_msg);
            }
        auto actual_range = offsetEnd()-offset();
        if(expected_range > actual_range)
            {
            auto error_msg = 
            format("SafePtr: expected_range > actual_range (%d > %d)",
                    expected_range,actual_range);
            throw std::runtime_error(error_msg);
            }
        return (p_+offset_);
        }

    explicit operator bool() const { return bool(p_); }

    SafePtr&
    operator+=(size_t shift)
        {
        //if(!p_) 
        //    {
        //    throw std::runtime_error("SafePtr: incrementing (+=) null pointer");
        //    }
        offset_ += shift;
        return *this;
        }

    SafePtr&
    operator-=(size_t shift)
        {
        //if(!p_) throw std::runtime_error("SafePtr: decrementing (-=) null pointer");
        offset_ -= shift;
        return *this;
        }

    SafePtr&
    operator++()
        {
        //if(!p_) throw std::runtime_error("SafePtr: incrementing null pointer");
        ++offset_;
        return *this;
        }

    SafePtr&
    operator++(int)
        {
        //if(!p_) throw std::runtime_error("SafePtr: incrementing null pointer");
        ++offset_;
        return *this;
        }

    SafePtr&
    operator--()
        {
        //if(!p_) throw std::runtime_error("SafePtr: decrementing null pointer");
        --offset_;
        return *this;
        }

    SafePtr&
    operator--(int)
        {
        //if(!p_) throw std::runtime_error("SafePtr: decrementing null pointer");
        --offset_;
        return *this;
        }

    reference
    operator*() { return *safeFront(); }

    pointer
    operator->() { return safeFront(); }

    reference
    operator[](size_t ind)
        {
        if(!p_) throw std::runtime_error("SafePtr operator[]: dereferencing null pointer");
        auto os = offset_+ind;
        if(os >= offset_end_)
            {
            auto error_msg = format("SafePtr operator[](ind=%d): (offset+ind) >= offset_end (%d >= %d)",ind,os,offset_end_);
            throw std::runtime_error(error_msg);
            }
        return *(p_+os);
        }

    bool
    operator!=(SafePtr const& other) const 
        { 
        if(p_ != other.p_)
            throw std::runtime_error("SafePtr: error, comparing two different starting pointers");
        return (offset_ != other.offset_);
        }

    //bool
    //operator!=(const_pointer other) const { return get() != other; }


    pointer
    getStart() const 
        { 
        return p_;
        }

    private:

    pointer
    safeFront() const 
        {
        if(!p_) throw std::runtime_error("SafePtr: dereferencing null pointer");
        if(!validOffset())
            {
            auto error_msg = 
            format("SafePtr: offset >= offset_end (%d >= %d)",
                   offset_,offset_end_);
            throw std::runtime_error(error_msg);
            }
        return (p_+offset_);
        }

    };

template<typename T>
SafePtr<T>
operator+(size_t inc, SafePtr<T> p) { p += inc; return p; }

template<typename T>
SafePtr<T>
operator+(SafePtr<T> p, size_t inc) { p += inc; return p; }

template<typename T>
SafePtr<T>
makeSafePtr(T* pt, size_t offset_end)
    {
    return SafePtr<T>(pt,offset_end);
    }

template<typename T>
SafePtr<T>
makeSafePtr(T* pt, size_t offset, size_t offset_end)
    {
    return SafePtr<T>(pt,offset,offset_end);
    }

template<typename NewType, typename OldType>
SafePtr<NewType>
reinterpret(SafePtr<OldType> const& p)
    {
    //if(!p.validOffset()) throw std::runtime_error("SafePtr: invalid offset found before attempting reinterpret_cast");
    auto nptr = reinterpret_cast<NewType*>(p.get());
    auto new_end = (p.range()*sizeof(OldType))/sizeof(NewType);
    return makeSafePtr(nptr,new_end);
    }


//use like SAFE_PTR_CHECK_SIZE(sp,assumed_size);
//checks if SafePtr sp has assumed_size elements in its range
#define SAFE_PTR_CHECK_SIZE(SP,SZ) assert(SP.validOffset() && ((SP.offsetEnd())-(SP.offset()))==SZ)

#ifdef DEBUG

//SafePtr versions of macros
#define MAKE_SAFE_PTR(P,SZ) makeSafePtr(P,SZ)
#define MAKE_SAFE_PTR_OFFSET(P,OFF,SZ) makeSafePtr(P,OFF,SZ)
#define SAFE_REINTERPRET(NT,SP) reinterpret<NT>(SP)
#define SAFE_PTR_GET(SP,SZ) SP.safeGet(SZ)
#define SAFE_PTR_OF(T) SafePtr<T>

#else

//bare pointer versions of macros
#define MAKE_SAFE_PTR(P,SZ) (P)
#define MAKE_SAFE_PTR_OFFSET(P,OFF,SZ) ((P)+(OFF))
#define SAFE_PTR_GET(P,SZ) P
#define SAFE_REINTERPRET(NT,SP) reinterpret_cast<NT*>(SP)
#define SAFE_PTR_OF(T) T*

#endif


} //namespace itensor

#endif
