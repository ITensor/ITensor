//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLETENSOR_H
#define __ITENSOR_SIMPLETENSOR_H

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
    T* data_ = nullptr; // always used to point to the storage
    storage vec_;      // used if it owns its storage
    public:

    // Constructor examples 1) simpletensor T{2,3,4}; 3 dimensional; first index from 0 to 1, etc
    //  if it owns storage  2) simpletensor T(v);  v a vector<int> or vector<long int>
    //                      3) simpletensor T;  uninitialized
    // resize has the same sets of arguments
    // To construct a simpletensor which does not own its data:
    //     simpletensor(T* dat, const vector<U>& dims)

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

    // U == int or long
    template<typename U, typename InputIterator>
    simpletensor(const std::vector<U>& dims,
                 InputIterator first,
                 InputIterator last)
        : 
        vec_(first,last)
        { 
        set_n(dims); 
        initNoAlloc(true); 
        }

    simpletensor(index&& dims,
                 storage&& vv)
        : 
        n_(std::move(dims)),
        vec_(std::move(vv))
        { 
        initNoAlloc(true); 
        }

    template<typename U> // U == int or long
    simpletensor(const T* dat, 
                 const std::vector<U>& dims) 
        : 
        data_(const_cast<T*>(dat))
        { 
        set_n(dims); 
        init(false); 
        }

    simpletensor(const T* dat, 
                 index&& dims) 
        : 
        n_(std::move(dims)),
        data_(const_cast<T*>(dat))
        { 
        init(false); 
        }

    // number of indices (tensor rank)
    long 
    r() const { return n_.size(); }

    long 
    size() const { return vec_.size(); }

    explicit operator bool() const { return bool(data_);}

    bool
    ownstorage() const { return (!vec_.empty() && data_==&vec_.front()); }

    // dimension of index i
    long 
    n(long i) const { return n_[i]; }

    // stride of index i
    long 
    stride(long i) const { return stride_[i]; }

    // direct access to element
    const T& 
    v(long i) const { return data_[i]; }

    // direct access to element
    T& 
    vref(long i) { return data_[i]; }

    // direct access to data
    const T*
    data() const { return data_; }

    // direct access to data
    T*
    data() { return data_; }

    // access storage
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

    template<typename Iterable>
    long
    indIterable(const Iterable& inds) const
        {
        long ii = 0, 
             i = 0;
        for(auto& j : inds)
            {
            ii += stride_[i] * j;
            ++i;
            }
        return ii;
        }

    template<typename U>
    long
    ind(const std::vector<U>& inds) const
        {
        return indIterable(inds);
        }

    template<typename U, size_t size>
    long
    ind(const std::array<U,size>& inds) const
        {
        return indIterable(inds);
        }

    template <typename... Inds>
    long
    ind(long i0, Inds... inds) { return __ind<0>(i0,inds...); }


    template <typename Ind0, typename... Inds>
    T& 
    operator()(const Ind0& i0, Inds... rest) { return data_[ind(i0,rest...)]; }

    template <typename Ind0, typename... Inds>
    const T& 
    operator()(Ind0 i0, Inds... rest) const { return data_[ind(i0,rest...)]; }

    void 
    clear()
        {
        vec_.clear();
        n_.clear();
        stride_.clear();
        data_ = nullptr;
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

    // assumes n_ has already been set
    long
    computeStride()
        {
        stride_.resize(r());
        long len = 1;
        for(long i = 0; i < r(); ++i)
            {
            stride_[i] = len;
            len *= n_[i];
            }
        return len;
        }

    // assumes n_ has already been set and
    // vec_ has already been allocated to correct size
    void 
    initNoAlloc(bool ownstore)
        {
        auto len = computeStride();
        if(ownstore && len > 0)
            {
            if(vec_.size() != len) throw std::runtime_error("Wrong size of input data");
            data_ = &vec_.front();
            }
        }

    // assumes n_ has already been set
    void 
    init(bool ownstore)
        {
        auto len = computeStride();
        if(ownstore && len > 0)
            {
            vec_.resize(len);
            data_ = &vec_.front();
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
