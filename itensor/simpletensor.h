//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLETENSOR_H
#define __ITENSOR_SIMPLETENSOR_H

#include "range.h"
#include "detail/gcounter.h"
#include <iostream>

namespace itensor {

template<typename T, typename RangeT = Range>
class tensorref
    {
    public:
    using value_type = T;
    using range_type = RangeT;
    using iterator = T*;
    using const_iterator = const T*;
    private:
    const range_type* inds_ = nullptr;  // size of each index;  first index dim is n[0]; indices from 0 to n[j]-1
    T* data_ = nullptr;                 // points to data storage
    public:

    tensorref() { }

    tensorref(const T* dat, 
              const RangeT& inds) 
        : 
        inds_(&inds),
        data_(const_cast<T*>(dat))
        { }

    void
    init(const T* dat,
         const RangeT& inds)
        {
        inds_ = &inds;
        data_ = const_cast<T*>(dat);
        }

    // number of indices (tensor rank)
    long 
    r() const { return inds_->size(); }

    long 
    size() const { return area(*inds_); }

    explicit operator bool() const { return bool(data_);}

    bool
    ownstorage() const { return false; }

    // dimension of index i
    long 
    n(long i) const { return (*inds_).dim(i); }

    // stride of index i
    long 
    stride(long i) const { return (*inds_).stride(i); }

    const range_type& 
    inds() const { return *inds_; }

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

    iterator
    begin() { return data_; }
    iterator
    end() { return data_+size(); }
    const_iterator
    begin() const { return data_; }
    const_iterator
    end() const { return data_+size(); }
    const_iterator
    cbegin() const { return begin(); }
    const_iterator
    cend() const { return end(); }

    template <typename... Inds>
    T& 
    operator()(const Inds&... ii) { return data_[ind(*inds_,ii...)]; }

    template <typename... Inds>
    const T& 
    operator()(const Inds&... ii) const { return data_[ind(*inds_,ii...)]; }

    void 
    clear()
        {
        inds_ = nullptr;
        data_ = nullptr;
        }

    void
    swap(tensorref& other)
        {
        std::swap(inds_,other.inds_);
        std::swap(data_,other.data_);
        }

    };

template<typename T, typename RangeT>
tensorref<T,RangeT>
make_tensorref(const T* dat,
               const RangeT& inds)
    {
    return tensorref<T,RangeT>(dat,inds); 
    }

template<typename T, typename RangeT = Range>
class tensor : public tensorref<T,RangeT>
    {
    public:
    using parent = tensorref<T,RangeT>;
    using value_type = typename parent::value_type;
    using storage_type = std::vector<T>;
    using range_type = typename parent::range_type;
    using iterator = typename parent::iterator;
    using const_iterator = typename parent::const_iterator;
    private:
    range_type inds_;   // size of each index;  first index dim is n[0]; indices from 0 to n[j]-1
    storage_type vec_;  // data storage
    public:

    tensor() : parent(nullptr,inds_) { }

    //constructor taking any number of index dimensions
    //e.g. tensor(4,4,3,4,2);
    template <typename... Dims>
    tensor(long d0, Dims... rest)
        :
        inds_(d0,rest...)
        {
        init();
        }

    tensor(std::initializer_list<long> dims)
        :
        inds_(dims)
        { 
        init(); 
        }

    template<typename U> // U == int or long
    tensor(const std::vector<U>& dims)
        :
        inds_(dims)
        { 
        init(); 
        }

    tensor(range_type&& dims)
        :
        inds_(std::move(dims))
        { 
        init(); 
        }

    tensor(range_type&& dims,
           storage_type&& v)
        :
        inds_(std::move(dims)),
        vec_(std::move(v))
        { 
        initNoAlloc(); 
        }

    // U == int or long
    template<typename U>
    tensor(const std::vector<U>& dims,
                 const T& val)
        : 
        inds_(dims)
        { 
        init(val);
        }

    // U == int or long
    template<typename U, typename InputIterator>
    tensor(const std::vector<U>& dims,
                 InputIterator first,
                 InputIterator last)
        : 
        inds_(dims),
        vec_(first,last)
        { 
        initNoAlloc();
        }

    template<typename InputIterator>
    tensor(const range_type& dims,
                 InputIterator first,
                 InputIterator last)
        : 
        inds_(dims),
        vec_(first,last)
        { 
        initNoAlloc();
        }

    template<typename InputIterator>
    tensor(range_type&& dims,
                 InputIterator first,
                 InputIterator last)
        : 
        inds_(std::move(dims)),
        vec_(first,last)
        { 
        initNoAlloc();
        }

    tensor(const tensor& other)
        {
        operator=(other);
        }

    tensor(tensor&& other)
        {
        operator=(std::move(other));
        }

    tensor&
    operator=(const tensor& other)
        { 
        inds_ = other.inds_;
        vec_ = other.vec_;
        if(other.ownstorage())
            {
            initNoAlloc();
            }
        else
            {
            throw std::runtime_error("Cannot currently assign tensorref to tensor");
            }
        return *this;
        }

    tensor&
    operator=(tensor&& other)
        { 
        inds_ = std::move(other.inds_);
        if(other.ownstorage())
            {
            vec_ = std::move(other.vec_);
            initNoAlloc();
            }
        else
            {
            throw std::runtime_error("Cannot currently assign tensorref to tensor");
            }
        return *this;
        }

    long 
    size() const { return vec_.size(); }

    bool
    ownstorage() const { return true; }

    // access storage
    const storage_type&
    store() const { return vec_; }
    storage_type&
    store() { return vec_; }

    // resize({2,4,3})
    void 
    resize(std::initializer_list<long> dims) 
        { 
        inds_ = range_type(dims);
        init(); 
        }

    template<typename U> 
    void 
    resize(const std::vector<U>& dims) 
        { 
        inds_ = range_type(dims);
        init(); 
        }

    void 
    clear()
        {
        inds_.clear();
        vec_.clear();
        parent::clear();
        }

    void
    swap(tensor& other)
        {
        inds_.swap(other.inds_);
        vec_.swap(other.vec_);
        parent::swap(other);
        }

    private:

    // assumes vec_ has already been allocated to correct size
    void 
    initNoAlloc()
        {
#ifdef DEBUG
        auto len = area(inds_);
        if(len == 0) throw std::runtime_error("Zero area in tensor");
        if(vec_.size() != size_t(len)) throw std::runtime_error("Wrong size of input data");
#endif
        parent::init(vec_.data(),inds_);
        }

    void 
    init()
        {
        auto len = area(inds_);
#ifdef DEBUG
        if(len == 0) throw std::runtime_error("Zero area in tensor");
#endif
        if(vec_.size() != size_t(len)) vec_ = storage_type(len);
        parent::init(vec_.data(),inds_);
        }

    void 
    init(const T& val)
        {
        auto len = area(inds_);
#ifdef DEBUG
        if(len == 0) throw std::runtime_error("Zero area in tensor");
#endif
        vec_ = storage_type(len,val);
        parent::init(vec_.data(),inds_);
        }
    };

template<typename T, typename RangeT, typename... Args>
long
ind(const tensorref<T,RangeT>& t,
    Args&&... args)
    {
    return ind(t.inds(),std::forward<Args>(args)...);
    }

template<typename T, typename RangeT, typename... Args>
long
ind(const tensor<T,RangeT>& t,
    Args&&... args)
    {
    return ind(t.inds(),std::forward<Args>(args)...);
    }

void
plusEqData(double fac,
           const double *d1,
           double *d2,
           int size);

// t2 += fac*t1, based on BLAS axpy
template<typename RangeT>
void
plusEq(double fac,
       const tensorref<double,RangeT>& t1,
       tensorref<double,RangeT>& t2)
    {
    plusEqData(fac,t1.data(),t2.data(),t1.size());
    }

namespace detail {
template<typename T, typename RangeT>
std::ostream&
printTensor(const char* type, std::ostream& s, const tensorref<T,RangeT>& t)
    {
    auto rank = t.r();
    auto gc = detail::GCounter(0,rank-1,0);
    s << type << "(";
    for(int i = 0; i < rank; ++i)
        {
        s << t.n(i);
        if(i != rank-1) s << ",";
        gc.setInd(i,0,t.n(i)-1);
        }
    s << ") = \n";
    for(; gc.notDone(); ++gc)
        {
        s << "  (";
        for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
            {
            s << +gc.i(ii);
            if(ii < gc.i.maxi()) s << ",";
            }
        s << ") ";

        s << t(gc.i) << "\n";
        }
    return s;
    }
}

template<typename T, typename RangeT>
std::ostream&
operator<<(std::ostream& s, const tensorref<T,RangeT>& t)
    {
    return detail::printTensor("tensorref",s,t);
    }

template<typename T, typename RangeT>
std::ostream&
operator<<(std::ostream& s, const tensor<T,RangeT>& t)
    {
    return detail::printTensor("tensor",s,t);
    }

} //namespace itensor

#endif
