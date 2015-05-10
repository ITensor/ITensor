//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SAFE_PTR_H
#define __ITENSOR_SAFE_PTR_H

#include "print.h"

namespace itensor {
namespace detail {

template<typename T>
class SafePtr
    {
    public:
    using value_type = std::remove_const_t<T>;
    using pointer = T*;
    using reference = T&;
    private:
    pointer p_ = nullptr;
    size_t offset_ = 0;
    size_t offset_end_ = 0;
    public:

    SafePtr(T* pt, size_t offset_end)
        : p_(pt), offset_(0), offset_end_(offset_end)
        { 
        if(!p_) throw std::runtime_error("SafePtr: pointer is null");
        }

    SafePtr(T* pt, size_t offset, size_t offset_end)
        : p_(pt), offset_(offset), offset_end_(offset_end)
        { 
        if(!p_) throw std::runtime_error("SafePtr: pointer is null");
        }

    explicit operator bool() const { return bool(p_); }

    SafePtr&
    operator+=(size_t shift)
        {
        if(!p_) throw std::runtime_error("SafePtr: incrementing null pointer");
        offset_ += shift;
        return *this;
        }

    SafePtr&
    operator++()
        {
        if(!p_) throw std::runtime_error("SafePtr: incrementing null pointer");
        ++offset_;
        return *this;
        }

    SafePtr&
    operator++(int)
        {
        if(!p_) throw std::runtime_error("SafePtr: incrementing null pointer");
        ++offset_;
        return *this;
        }

    reference
    operator*()
        {
        if(!p_) throw std::runtime_error("SafePtr: dereferencing null pointer");
        if(offset_ >= offset_end_) 
            {
            auto error_msg = format("SafePtr: offset >= offset_end (%d >= %d)",offset_,offset_end_);
            throw std::runtime_error(error_msg);
            }
        return *(p_+offset_);
        }
    };

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


} //namespace detail
} //namespace itensor

#endif
