//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TENSOR_TYPES_H
#define __ITENSOR_TENSOR_TYPES_H

#include <stdexcept>
#include <iostream>
#include "itensor/util/infarray.h"
#include "itensor/util/vararray.h"
#include "itensor/util/timers.h"

namespace itensor {

using Label = InfArray<long,11ul>; //sizeof(InfArray<long,11ul>)==128
//using Label = VarArray<long,15ul>; //sizeof(VarArray<long,15ul>)==128
//using Label = VarArray<long,31ul>; //sizeof(VarArray<long,31ul>)==256

// Storage for TenRef
template<typename T>
class DataRange;

using Data  = DataRange<double>;
using cData = DataRange<const double>;

template<typename T>
class DataRange
    {
    public:
    using value_type = typename std::decay<T>::type;
    using pointer = T*;
    using reference = T&;
    using const_reference = const T&;
    private:
    pointer pdata_;
    size_t  size_;
    public:

    DataRange() : pdata_(nullptr), size_(0) { }

    DataRange(pointer pdat,
              size_t  size)
      : pdata_(pdat),
        size_(size)
        { }

    DataRange(pointer pdat,
              size_t  offset,
              size_t  size)
      : pdata_(pdat+offset),
        size_(size-offset)
        { }

    explicit
    operator bool() const { return static_cast<bool>(pdata_); }

    operator DataRange<const T>() const 
        {
        return DataRange<const T>{pdata_,size_};
        }

    reference
    operator[](size_t n)
        {
#ifdef DEBUG
        if(n >= size_) 
            {
            std::cout << "offset = " << n << std::endl;
            std::cout << "max offset (size-1) = " << (size_==0 ? -1ul : size_-1ul) << std::endl;
            Error("offset exceeded size in DataRange operator[]");
            }
#endif
        return pdata_[n];
        }

    const_reference
    operator[](size_t n) const
        {
#ifdef DEBUG
        if(n >= size_) 
            {
            std::cout << "offset = " << n << std::endl;
            std::cout << "max offset (size-1) = " << (size_==0 ? -1ul : size_-1ul) << std::endl;
            Error("offset exceeded size in DataRange operator[]");
            }
#endif
        return pdata_[n];
        }


    pointer
    data() const { return pdata_; }

    size_t
    size() const { return size_; }

    DataRange
    operator+(size_t off) const
        {
        auto res = DataRange{};
        res.pdata_ = pdata_+off;
#ifdef DEBUG
        if(off > size_)
            {
            Error("attempt to add offset to data greater than size");
            }
#endif
        res.size_ = size_-off;
        return res;
        }

    void
    clear()
        {
        pdata_ = nullptr;
        size_  = 0;
        }

    DataRange<value_type>
    cast_away_const() const
        {
        return DataRange<value_type>{const_cast<value_type*>(pdata_),size_};
        }
    };


template<typename T, size_t N>
std::ostream& 
operator<<(std::ostream & s, InfArray<T,N> const& v)
    {
    if(v.empty()) return s;
    decltype(v.size()) j = 0;
    for(; 1+j < v.size(); ++j)
        {
        s << v[j] << ",";
        }
    s << v[j];
    return s;
    }

//template<typename T>
//std::ostream& 
//operator<<(std::ostream & s, std::vector<T> const& v)
//    {
//    if(v.empty()) return s;
//    decltype(v.size()) j = 0;
//    for(; 1+j < v.size(); ++j)
//        {
//        s << v[j] << ",";
//        }
//    s << v[j];
//    return s;
//    }

} //namespace itensor

#endif
