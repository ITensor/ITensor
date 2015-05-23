//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEN_H_
#define __ITENSOR_TEN_H_

#include "itensor/util/types.h"
#include "itensor/tensor/range.h"

namespace itensor {

template<typename T, typename RangeT>
class Tensor;

template<typename T, typename RangeT>
class TensorRef;

using Ten     = Tensor<Real,Range>;
template<typename RangeT>
using TenRef  = TensorRef<Real,RangeT>;
template<typename RangeT>
using TenRefc = TensorRef<const Real,RangeT>;

template<typename T, typename RangeT = Range>
class TensorRef
    { 
    public:
    using iterator = T*;
    using const_iterator = const T*;
    using value_type = std::remove_const_t<T>;
    using pointer = T*;
    using reference = T&;
    using size_type = long;
    using range_type = RangeT;
    using tensor_type = std::conditional_t<std::is_const<T>::value,
                                           const Tensor<value_type,range_type>,
                                           Tensor<value_type,range_type>>;
    private:
    pointer pdata_ = nullptr;
    const range_type* prange_ = nullptr;
    public:

    TensorRef() { }

    TensorRef(pointer pdat,
              const range_type& range)
        :
        pdata_(pdat),
        prange_(&range)
        { }

    TensorRef(tensor_type& t) { pointTo(t); }

    TensorRef&
    operator=(tensor_type& t) { pointTo(t); return *this; }

    TensorRef(tensor_type&& t) = delete;

    TensorRef&
    operator=(tensor_type&& t) = delete;

    operator TensorRef<const value_type,range_type>() const { return TensorRef<const value_type,range_type>(data(),range()); }

    size_type
    r() const { return prange_->r(); }

    size_type 
    size() const { return area(*prange_); }

    explicit operator bool() const { return bool(pdata_);}

    size_type
    dim(size_type i) const { return prange_->dim(i); }

    size_type
    stride(size_type i) const { return prange_->stride(i); }

    const range_type&
    range() const { return *prange_; }

    // direct access to data
    pointer
    data() const { return pdata_; }

    template <typename... Inds>
    reference
    operator()(Inds&&... ii) const { return pdata_[ind(*prange_,std::forward<Inds>(ii)...)]; }

    void
    clear() { pdata_ = nullptr; prange_ = nullptr; }

    //Currently assuming TensorRefs point to contiguous data:

    iterator
    begin() const { return pdata_; }

    iterator
    end() const { return pdata_+size(); }

    const_iterator
    cbegin() const { return pdata_; }

    const_iterator
    cend() const { return pdata_+size(); }

    private:
    void
    pointTo(tensor_type& t);
    };

template<typename T, typename RangeT>
auto
makeTensorRef(T* p,
              const RangeT& range)
    {
    return TensorRef<T,RangeT>(p,range);
    }

template<typename T, typename RangeT = Range>
class Tensor
    {
    public:
    using storage_type = std::vector<T>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;
    using value_type = std::remove_const_t<T>;
    using pointer = std::add_pointer_t<value_type>;
    using const_pointer = std::add_pointer_t<const value_type>;
    using reference = std::add_lvalue_reference_t<value_type>;
    using const_reference = const reference;
    using size_type = long;
    using range_type = RangeT;
    public:
    range_type range_;
    storage_type data_;
    public:

    Tensor() { }

    template<typename... Dims>
    Tensor(size_type d0, Dims&&... rest)
        : range_(d0,std::forward<Dims>(rest)...)
        {
        init();
        }

    Tensor(const range_type& range)
        : range_(range)
        {
        init();
        }

    Tensor(range_type&& range)
        : range_(std::move(range))
        {
        init();
        }

    Tensor(range_type&& range,
           storage_type&& data)
        : 
        range_(std::move(range)),
        data_(std::move(data))
        { 
#ifdef DEBUG
        auto len = area(range_);
        if(len==0) throw std::runtime_error("Zero are in Tensor");
        if(data_.size() != size_t(len)) throw std::runtime_error("Wrong size of input data");
#endif
        }

    Tensor(const Tensor& other) { assignFrom(other); }

    Tensor(Tensor&& other) { moveFrom(std::move(other)); }

    Tensor&
    operator=(const Tensor& other) { assignFrom(other); return *this; }

    Tensor& 
    operator=(Tensor&& other) { moveFrom(std::move(other)); return *this; }

    explicit
    Tensor(const TensorRef<const value_type,range_type>& ref) { assignFromRef(ref); }

    explicit
    Tensor(const TensorRef<value_type,range_type>& ref) { assignFromRef(ref); }

    Tensor&
    operator=(const TensorRef<const value_type>& ref) { assignFromRef(ref); return *this; }

    Tensor&
    operator=(const TensorRef<value_type>& ref) { assignFromRef(ref); return *this; }

    explicit operator bool() const { return !data_.empty(); }

    size_type
    r() const { return range_.r(); }

    size_type
    size() const { return data_.size(); }

    size_type
    dim(size_type i) const { return range_.dim(i); }

    size_type
    stride(size_type i) const { return range_.stride(i); }

    const range_type&
    range() const { return range_; }

    template <typename... Inds>
    reference
    operator()(Inds&&... ii) { return data_[ind(range_,std::forward<Inds>(ii)...)]; }

    template <typename... Inds>
    const_reference
    operator()(Inds&&... ii) const { return data_[ind(range_,std::forward<Inds>(ii)...)]; }

    iterator
    begin() { return data_.begin(); }

    iterator
    end() { return data_.end(); }

    const_iterator
    begin() const { return data_.cbegin(); }

    const_iterator
    end() const { return data_.cend(); }

    const_iterator
    cbegin() const { return data_.cbegin(); }

    const_iterator
    cend() const { return data_.cend(); }

    pointer
    data() { return data_.data(); }

    const_pointer
    data() const { return data_.data(); }

    void
    clear() { data_.clear(); range_.clear(); }

    storage_type&
    store() { return data_; }

    const storage_type&
    store() const { return data_; }

    private:

    void
    init()
        {
        auto len = area(range_);
#ifdef DEBUG
        if(len == 0) throw std::runtime_error("Zero area in tensor");
#endif
        data_.assign(len,0.);
        }

    void
    assignFromRef(const TensorRef<const value_type>& other)
        {
        //Copy data from other contiguously into data_
        data_ = storage_type(other.cbegin(),other.cend());
        }

    void
    assignFrom(const Tensor& other)
        {
        if(&other == this) return;
        range_ = other.range_;
        data_ = other.data_;
        }

    void
    moveFrom(Tensor&& other)
        {
        range_ = std::move(other.range_);
        data_ = std::move(other.data_);
        other.clear();
        }
    };

template<typename T, typename R>
void TensorRef<T,R>::
pointTo(tensor_type& t)
    {
    pdata_ = t.data();
    prange_ = &t.range();
    }


//
// makeRef functions
//

template<typename T, typename R>
auto
makeRef(TensorRef<T,R>& t) { return t; }

template<typename T, typename R>
auto
makeRef(TensorRef<const T,R>& t) { return t; }

template<typename T, typename R>
auto
makeRef(Tensor<T,R>& t) { return TensorRef<T,R>(t); }

template<typename T, typename R>
auto
makeRef(const Tensor<T,R>& t) { return TensorRef<const T>(t); }

//This version of makeRef intended to fail,
//forbids explicitly making TensorRefs to temporaries
template<typename T, typename R, typename... Rest>
auto
makeRef(Tensor<T,R>&& t, Rest&&... args) { return TensorRef<T,R>(std::move(t)); }

template<typename T, typename R>
auto
makeRefc(TensorRef<T,R>& t) { return TensorRef<const T,R>(t); }

template<typename T, typename R>
auto
makeRefc(TensorRef<const T,R>& t) { return t; }

template<typename T, typename R>
auto
makeRefc(Tensor<T,R>& t) { return TensorRef<const T,R>(t); }

template<typename T, typename R>
auto
makeRefc(const Tensor<T,R>& t) { return TensorRef<const T>(t); }

//This version of makeRefc intended to fail,
//forbids explicitly making TensorRefs to temporaries
template<typename T, typename R, typename... Rest>
auto
makeRefc(Tensor<T,R>&& t, Rest&&... args) { return TensorRef<const T,R>(std::move(t)); }

} //namespace itensor

#endif
