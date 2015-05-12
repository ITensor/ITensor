//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SAFE_PTR_H
#define __ITENSOR_SAFE_PTR_H

#include "print.h"

namespace itensor {

template<typename T>
class SafePtr
    {
    public:
    using value_type = std::remove_const_t<T>;
    using pointer = T*;
    using reference = T&;
    using const_pointer = const T*;
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

    size_t
    offset() const { return offset_; }
    size_t
    offsetEnd() const { return offset_end_; }
    bool
    validOffset() const { return (offset_ < offset_end_); }

    pointer
    get() const 
        { 
        if(!p_) return nullptr;
        return p_+offset_; 
        }

    pointer
    safeGet() const 
        {
        if(!p_) throw std::runtime_error("SafePtr: dereferencing null pointer");
        if(!validOffset())
            {
            auto error_msg = format("SafePtr: offset >= offset_end (%d >= %d)",offset_,offset_end_);
            throw std::runtime_error(error_msg);
            }
        return (p_+offset_);
        }

    explicit operator bool() const { return bool(p_); }

    SafePtr&
    operator+=(size_t shift)
        {
        if(!p_) throw std::runtime_error("SafePtr: incrementing (+=) null pointer");
        offset_ += shift;
        return *this;
        }

    SafePtr&
    operator-=(size_t shift)
        {
        if(!p_) throw std::runtime_error("SafePtr: decrementing (-=) null pointer");
        offset_ -= shift;
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

    SafePtr&
    operator--()
        {
        if(!p_) throw std::runtime_error("SafePtr: decrementing null pointer");
        --offset_;
        return *this;
        }

    SafePtr&
    operator--(int)
        {
        if(!p_) throw std::runtime_error("SafePtr: decrementing null pointer");
        --offset_;
        return *this;
        }

    reference
    operator*() { return *safeGet(); }

    pointer
    operator->() { return safeGet(); }

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
    operator!=(const SafePtr& other) const { return get() != other.get(); }

    bool
    operator!=(const_pointer other) const { return get() != other; }

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

#ifdef DEBUG
#define MAKE_SAFE_PTR(X,Y) makeSafePtr(X,Y)
#define MAKE_SAFE_PTR3(X,Y,Z) makeSafePtr(X,Y,Z)
#else
#define MAKE_SAFE_PTR(X,Y) (X)
#define MAKE_SAFE_PTR3(X,Y,Z) ((X)+(Y))
#endif


} //namespace itensor

#endif
