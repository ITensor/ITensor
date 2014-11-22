//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLETENSOR_H
#define __ITENSOR_SIMPLETENSOR_H

#include "matrixref.h"
#include "autovector.h"

namespace itensor {

//
// simpletensor SRW 8/25/14
//
// simple tensor class for any type of element

// Possible optimizations/redesigns:
// o Don't store stride_, instead compute on the fly
//   (may want to keep stride_ for future slicing operations, 
//   or create a less-simple tensor for that)

// TODO
//
// o Move ind logic either to Range or as external function
// 

class Range
    {
    public:
    struct index 
        { 
        long dim = 0, 
             stride = 0; 
        index() { }
        index(long dim_, long stride_) : dim(dim_),stride(stride_) { }
        explicit operator long() const { return dim; }
        };
    using storage = std::vector<index>;

    Range() { }

    template <typename... Dims>
    Range(long d0, Dims... rest)
        {
        const auto ndim = 1+sizeof...(rest);
        std::array<long,ndim> dims = {{d0,static_cast<long>(rest)...}};
        init(dims);
        }

    template <typename Iterable>
    Range(const Iterable& v)
        {
        init(v);
        }

    const index&
    operator[](int i) const { return inds_[i]; }
    index&
    operator[](int i) { return inds_[i]; }

    const index&
    at(int i) const { return inds_.at(i); }
    index&
    at(int i) { return inds_.at(i); }

    size_t
    size() const { return inds_.size(); }

    void
    clear() { return inds_.clear(); }

    bool
    empty() const { return inds_.empty(); }

    index&
    front() { return inds_.front(); }
    const index&
    front() const { return inds_.front(); }

    index&
    back() { return inds_.back(); }
    const index&
    back() const { return inds_.back(); }

    private:

    template<typename Iterable>
    void 
    init(const Iterable& v)
        {
        inds_.resize(v.size());
        long len = 1;
        for(size_t i = 0; i < v.size(); ++i)
            {
            inds_[i] = index(long{v[i]},len);
            len *= inds_[i].dim;
            }
        }

    storage inds_;
    };

long inline
area(const Range& r)
    { 
    if(r.empty()) return 0;
    return r.back().dim*r.back().stride;
    }


template<typename T>
class simpletensor
    {
    public:
    using value_type = T;
    using storage = std::vector<T>;
    using index_type = Range;
    private:
    index_type inds_;          // size of each index;  first index dim is n[0]; indices from 0 to n[j]-1
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
        inds_{d0,rest...}
        {
        init(true);
        }

    simpletensor(std::initializer_list<long> dims)
        :
        inds_(dims)
        { 
        init(true); 
        }

    template<typename U> // U == int or long
    simpletensor(const std::vector<U>& dims)
        :
        inds_(dims)
        { 
        init(true); 
        }

    simpletensor(Range&& dims)
        :
        inds_(std::forward<Range>(dims))
        { 
        init(true); 
        }

    // U == int or long
    template<typename U, typename InputIterator>
    simpletensor(const std::vector<U>& dims,
                 InputIterator first,
                 InputIterator last)
        : 
        inds_(dims),
        vec_(first,last)
        { 
        initNoAlloc(true);
        }

    template<typename InputIterator>
    simpletensor(Range&& dims,
                 InputIterator first,
                 InputIterator last)
        : 
        inds_(std::move(dims)),
        vec_(first,last)
        { 
        initNoAlloc(true);
        }

    template<typename U> // U == int or long
    simpletensor(const T* dat, 
                 const std::vector<U>& dims) 
        : 
        inds_(dims),
        data_(const_cast<T*>(dat))
        { 
        init(false); 
        }

    simpletensor(const T* dat, 
                 Range&& dims) 
        : 
        inds_(std::move(dims)),
        data_(const_cast<T*>(dat))
        { 
        init(false); 
        }

    // number of indices (tensor rank)
    long 
    r() const { return inds_.size(); }

    long 
    size() const { return area(inds_); }

    explicit operator bool() const { return bool(data_);}

    bool
    ownstorage() const { return (!vec_.empty() && data_==&vec_.front()); }

    // dimension of index i
    long 
    n(long i) const { return inds_[i].dim; }

    // stride of index i
    long 
    stride(long i) const { return inds_[i].stride; }

    // stride of index i
    index_type 
    inds() const { return inds_; }

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

    T*
    begin() { return data_; }
    T*
    end() { return data_+size(); }
    const T*
    begin() const { return data_; }
    const T*
    end() const { return data_+size(); }
    const T*
    cbegin() const { return begin(); }
    const T*
    cend() const { return end(); }

    VectorRefNoLink 
    vectorref() { return VectorRefNoLink(data,size()); }

    // resize({2,4,3})
    void 
    resize(std::initializer_list<long> dims) 
        { 
        init(true,dims); 
        }

    template<typename U> 
    void 
    resize(const std::vector<U>& dims) 
        { 
        init(true,dims); 
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

    template<typename U>
    long
    ind(const autovector<U>& inds) const
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
    operator()(const Ind0& i0, Inds... rest) const { return data_[ind(i0,rest...)]; }

    void 
    clear()
        {
        vec_.clear();
        inds_.clear();
        data_ = nullptr;
        }

    private:

    template<typename Iterable>
    long
    indIterable(const Iterable& inds) const
        {
        long ii = 0, 
             i = 0;
        for(auto& j : inds)
            {
            ii += inds_[i].stride * j;
            ++i;
            }
        return ii;
        }


    // assumes vec_ has already been allocated to correct size
    void 
    initNoAlloc(bool ownstore)
        {
        auto len = area(inds_);
        if(ownstore && len > 0)
            {
            if(vec_.size() != len) throw std::runtime_error("Wrong size of input data");
            data_ = &vec_.front();
            }
        }

    void 
    init(bool ownstore, T val = 0)
        {
        auto len = area(inds_);
        if(ownstore && len > 0)
            {
            vec_ = storage(len,val);
            data_ = &vec_.front();
            }
        }

    template <long i, typename... Inds>
    long
    __ind(long first, Inds... rest) const
        {
        return first*inds_[i].stride + __ind<i+1>(rest...);
        }

    template <long i>
    long
    __ind(long first) const
        {
        return first*inds_[i].stride;
        }

    };


}; //namespace itensor

#endif
