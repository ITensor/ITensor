//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLETENSOR_H
#define __ITENSOR_SIMPLETENSOR_H

#include <array>
#include "matrixref.h"

namespace itensor {

//
// simpletensor SRW 8/25/14
//
// simple tensor class for any type of element

// Possible optimizations/redesigns:
// o Don't store stride_, instead compute on the fly
//   (may want to keep stride_ for future slicing operations, 
//   or create a less-simple tensor for that)

template<typename T>
class simpletensor
    {
    public:
    using value_type = T;
    using storage = std::vector<T>;
    using index = std::vector<long>;
    private:
    index n_;          // size of each index;  first index dim is n[0]; indices from 0 to n[j]-1
    index stride_;     // spacing in memory for that index, stride_[0]==1 (column-major order)
    T* data = nullptr; // always used to point to the storage
    storage vec_;      // used if it owns its storage
    public:

    // Constructor examples 1) simpletensor T{2,3,4}; 3 dimensional; first index from 0 to 1, etc
    //  if it owns storage  2) simpletensor T(v);  v a vector<int> or vector<long int>
    //                      3) simpletensor T;  uninitialized
    // resize has the same sets of arguments
    // The only way to construct it where it does not own its data
    // is  simpletensor(T* dat, const vector<U>& dims)

    simpletensor() { }

    //constructor taking any number of index dimensions
    //e.g. simpletensor(4,4,3,4,2);
    template <typename... Dims>
    simpletensor(long d0, Dims... rest)
        :
        n_{d0,rest...}
        {
        init(true);
        }

    simpletensor(std::initializer_list<long> dims)
        { 
        n_ = dims; 
        init(true); 
        }

    template<typename U> // U == int or long
    simpletensor(const std::vector<U>& dims)
        { 
        set_n(dims); 
        init(true); 
        }

    template<typename U> // U == int or long
    simpletensor(const T* dat, 
                 const std::vector<U>& dims) 
        : 
        data(const_cast<T*>(dat))
        { 
        set_n(dims); 
        init(false); 
        }

    simpletensor(const T* dat, index&& dims) 
        : 
        n_(std::move(dims)),
        data(const_cast<T*>(dat))
        { 
        init(false); 
        }

    // number of indices (tensor rank)
    long 
    r() const { return n_.size(); }

    // dimension of index i
    long 
    n(long i) const { return n_[i]; }

    // stride of index i
    long 
    stride(long i) const { return stride_[i]; }

    long 
    size() const { return vec_.size(); }

    // direct access to data
    const T& 
    v(long i) const { return data[i]; }

    // direct access to data
    T& 
    vref(long i) { return data[i]; }

    const storage&
    store() const { return vec_; }

    VectorRefNoLink 
    vectorref() { return VectorRefNoLink(data,size()); }

    // resize({2,4,3})
    void 
    resize(std::initializer_list<long> dims) 
        { 
        n_ = dims; 
        init(true); 
        }

    template<typename U> 
    void 
    resize(const std::vector<U>& dims) 
        { 
        set_n(dims); 
        init(true); 
        }

    template<typename U>
    long 
    ind(const std::vector<U>& inds) const
        {
        long ii = 0, i = 0;
        for(auto& j : inds)
            {
            ii += stride_[i] * j;
            ++i;
            }
        return ii;
        }

    long 
    ind(std::initializer_list<long> inds) const
        {
        long ii = 0, i = 0;
        for(auto& j : inds)
            {
            ii += stride_[i] * j;
            ++i;
            }
        return ii;
        }

    long 
    ind(long i0) const { return i0; }

    template <typename... Inds>
    long
    ind(long i0, Inds... inds) { return __ind<0>(i0,inds...); }

    // Element access through t({i,j,k,l}) etc. any number of indices
    //                through t(i,j,k,l) etc. probably fastest
    //                through t(vector<long>& vi) generally
    T& 
    operator()(std::initializer_list<long> inds) { return data[ind(inds)]; }

    const T& 
    operator()(const index& inds) const { return data[ind(inds)]; }

    T& 
    operator()(long i0) { return data[ind(i0)]; }

    const T& 
    operator()(long i0) const { return data[ind(i0)]; }

    template <typename... Inds>
    T& 
    operator()(long i0, Inds... rest) { return data[ind(i0,rest...)]; }

    template <typename... Inds>
    const T& 
    operator()(long i0, Inds... rest) const { return data[ind(i0,rest...)]; }

    bool
    ownstorage() const { return (!vec_.empty() && data==&vec_.front()); }

    void 
    clear()
        {
        vec_.clear();
        n_.clear();
        stride_.clear();
        data = nullptr;
        }

    private:

    template<typename U>
    void 
    set_n(const std::vector<U>& v)
        {
        n_.resize(v.size());
        for(int i = 0; i < v.size(); ++i)
            n_[i] = long{v[i]};
        }

    void 
    init(bool ownstore) // if member-data n_ is initialized, take care of stride and data
        {
        stride_.resize(r());
        long len = 1;
        for(long i = 0; i < r(); ++i)
            {
            stride_[i] = len;
            len *= n_[i];
            }
        if(ownstore && len > 0)
            {
            vec_.resize(len);
            data = &vec_.front();
            }
        }

    template <long i, typename... Inds>
    long
    __ind(long first, Inds... rest)
        {
        return first*stride_[i] + __ind<i+1>(rest...);
        }

    template <long i>
    long
    __ind(long first)
        {
        return first*stride_[i];
        }

    };

}; //namespace itensor
#endif
