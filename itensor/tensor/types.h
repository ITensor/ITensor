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
#include "itensor/types.h"

namespace itensor {

template<typename T>
bool inline constexpr
isReal() { return std::is_same<stdx::decay_t<T>,Real>::value; }

template<typename T>
bool inline constexpr
isCplx() { return std::is_same<stdx::decay_t<T>,Cplx>::value; }

using IntArray = InfArray<long,11ul>; //sizeof(InfArray<long,11ul>)==128
using Labels = IntArray;
//using Labels = VarArray<long,15ul>; //sizeof(VarArray<long,15ul>)==128
//using Labels = VarArray<long,31ul>; //sizeof(VarArray<long,31ul>)==256

// Storage for TenRef
template<typename T>
class DataRange;

using Data  = DataRange<Real>;
using Datac = DataRange<const Real>;

using CData = DataRange<Cplx>;
using CDatac = DataRange<const Cplx>;

template<typename T>
class DataRange
    {
    public:
    using value_type = typename std::remove_reference<T>::type;
    using pointer = value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;
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
            if(size_ == 0) std::cout << "size_ = " << size_ << std::endl;
            else std::cout << "max offset (size-1) = " << (size_-1ul) << std::endl;
            Error("tensor (or vector/matrix) out of bounds access");
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
            if(size_ == 0) std::cout << "size_ = " << size_ << std::endl;
            else std::cout << "max offset (size-1) = " << (size_-1ul) << std::endl;
            Error("tensor (or vector/matrix) out of bounds access");
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

template<typename T>
DataRange<T>
makeDataRange(T * p, size_t size)
    {
    return DataRange<T>(p,size);
    }

template<typename T>
DataRange<const T>
makeDataRange(T const* p, size_t size)
    {
    return DataRange<const T>(p,size);
    }

template<typename T>
DataRange<T>
makeDataRange(T * p, size_t offset, size_t size)
    {
    return DataRange<T>(p,offset,size);
    }

template<typename T>
DataRange<const T>
makeDataRange(T const* p, size_t offset, size_t size)
    {
    return DataRange<const T>(p,offset,size);
    }

template<typename T>
DataRange<T>
sliceData(DataRange<T> d, size_t begin, size_t end)
    {
    auto size = end-begin;
#ifdef DEBUG
    if(begin > end) Error("begin > end in sliceData");
    if(begin+size > d.size()) 
        {
        println("d.size() = ",d.size());
        println("begin = ",begin);
        println("end = ",end);
        Error("sliceData invalid begin or end");
        }
#endif
    auto pb = d.data()+begin;
    return DataRange<T>(pb,size);
    }


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
