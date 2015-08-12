//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEN_H_
#define __ITENSOR_TEN_H_

#include "itensor/util/stdx.h"
#include "itensor/tensor/teniter.h"
#include "itensor/tensor/range.h"

namespace itensor {

template<typename T, typename RangeT>
class Ten;
template<typename T, typename RangeT>
class TenRef;

//Specialize to Real elements
template<typename RangeT>
using RTen     = Ten<Real,RangeT>;
template<typename RangeT>
using RTenRef  = TenRef<Real,RangeT>;
template<typename RangeT>
using RTenRefc = TenRef<const Real,RangeT>;

//Specialize to Real elements, RangeT=Range
using Tensor     = RTen<Range>;
using TensorRef  = RTenRef<Range>;
using TensorRefc = RTenRefc<Range>;

//Specialize to Cplx elements
//using CTen     = Ten<Cplx,RangeT>;
//template<typename RangeT>
//using CTenRef  = TenRef<Cplx,RangeT>;
//template<typename RangeT>
//using CTenRefc = TenRef<const Cplx,RangeT>;

//Specialize to Cplx elements, RangeT=Range
//using CTensor     = Tensor<Cplx,Range>;
//using CTensorRef  = TenRef<Cplx,Range>;
//using CTensorRefc  = TenRef<const Cplx,Range>;

template<typename T, typename RangeT>
class TenRef
    { 
    public:
    using value_type = std::decay_t<T>;
    using iterator = TenIter<T*,RangeT>;
    using const_iterator = TenIter<const T*,RangeT>;
    using pointer = T*;
    using reference = T&;
    using size_type = long;
    using range_type = std::remove_const_t<RangeT>;
    using tensor_type = std::conditional_t<std::is_const<T>::value,
                                           const Ten<value_type,range_type>,
                                           Ten<value_type,range_type>>;
    private:
    pointer pdata_ = nullptr;
    const range_type* prange_ = nullptr;
    range_type range_;
    public:

    TenRef() { }

    TenRef(pointer pdat,
           range_type const& range)
      : pdata_(pdat),
        prange_(&range)
        { }

    TenRef(pointer pdat,
           range_type && range)
      : pdata_(pdat),
        range_(std::move(range))
        { 
        prange_ = &range_;
        }

    TenRef(pointer pdat,
           const range_type* prange)
      : pdata_(pdat),
        prange_(prange)
        { }

    TenRef(TenRef const& t)
        {
        operator=(t);
        }

    TenRef&
    operator=(TenRef const& t);

    TenRef(TenRef && t)
        {
        operator=(std::move(t));
        }

    TenRef&
    operator=(TenRef && t);

    TenRef(tensor_type& t) { pointTo(t); }

    TenRef&
    operator=(tensor_type& t) { pointTo(t); return *this; }

    TenRef(tensor_type&& t) = delete;

    TenRef&
    operator=(tensor_type&& t) = delete;

    operator TenRef<const value_type,range_type>() const;

    bool
    ownRange() const { return prange_ == &range_; }

    size_type
    r() const { return prange_->r(); }

    size_type 
    size() const { return area(*prange_); }

    explicit operator bool() const { return bool(pdata_);}

    size_type
    extent(size_type i) const { return prange_->extent(i); }

    size_type
    stride(size_type i) const { return prange_->stride(i); }

    range_type const&
    range() const { return *prange_; }

    // direct access to data
    pointer
    data() const { return pdata_; }
    
    reference
    operator()() const;

    template <typename... Inds>
    reference
    operator()(Inds&&... ii) const;

    template<typename Indices>
    reference
    operator()(Indices const& ii) const;

    void
    clear() { pdata_ = nullptr; prange_ = nullptr; }

    iterator
    begin() const { return iterator(pdata_,*prange_); }

    iterator
    end() const { return iterator::makeEnd(*prange_); }

    const_iterator
    cbegin() const { return const_iterator(pdata_,*prange_); }

    const_iterator
    cend() const { return const_iterator::makeEnd(*prange_); }

    private:
    void
    pointTo(tensor_type& t);
    };

//Assign to referenced data
void 
operator&=(TensorRef const& a, TensorRefc b);

//Assign to referenced data
void
operator&=(TensorRef const& a, Tensor const& t);

template<typename T, typename RangeT>
auto
makeTenRef(T* p,
           const RangeT* prange)
    {
    using R = std::decay_t<std::remove_pointer_t<RangeT>>;
    //static_assert(!std::is_pointer<R>::value,"Error: RangeT is of pointer type");
    return TenRef<T,R>(p,prange);
    }

template<typename T, typename RangeT,
         class = std::enable_if_t<std::is_rvalue_reference<RangeT&&>::value
                               && !std::is_pointer<RangeT>::value> >
auto
makeTenRef(T* p,
           RangeT && range)
    {
    using R = std::decay_t<RangeT>;
    static_assert(!std::is_pointer<R>::value,"Error: RangeT is of pointer type");
    return TenRef<T,R>(p,std::move(range));
    }

template<typename T, typename RangeT>
class Ten
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

    Ten() { }

    template<typename... Dims>
    Ten(size_type d0, Dims&&... rest)
      : range_(d0,std::forward<Dims>(rest)...)
        {
        init();
        }

    Ten(range_type const& range)
      : range_(range)
        {
        init();
        }

    Ten(range_type && range)
      : range_(std::move(range))
        {
        init();
        }

    Ten(range_type && range,
        storage_type && data)
      : range_(std::move(range)),
        data_(std::move(data))
        { 
#ifdef DEBUG
        auto len = area(range_);
        if(len==0) throw std::runtime_error("Zero area in Ten");
        if(data_.size() != size_t(len)) throw std::runtime_error("Wrong size of input data");
#endif
        }

    Ten(Ten const& other) { assignFrom(other); }

    Ten(Ten && other) { moveFrom(std::move(other)); }

    Ten&
    operator=(Ten const& other) { assignFrom(other); return *this; }

    Ten& 
    operator=(Ten && other) { moveFrom(std::move(other)); return *this; }

    template<typename R>
    explicit
    Ten(TenRef<const value_type,R> const& ref) { assignFromRef(ref); }

    template<typename R>
    explicit
    Ten(TenRef<value_type,R> const& ref) { assignFromRef(ref); }

    template<typename R>
    Ten&
    operator=(TenRef<const value_type,R> const& ref) { assignFromRef(ref); return *this; }

    template<typename R>
    Ten&
    operator=(TenRef<value_type,R> const& ref) { assignFromRef(ref); return *this; }

    explicit operator bool() const { return !data_.empty(); }

    size_type
    r() const { return range_.r(); }

    size_type
    size() const { return data_.size(); }

    size_type
    extent(size_type i) const { return range_.extent(i); }

    size_type
    stride(size_type i) const { return range_.stride(i); }

    range_type const&
    range() const { return range_; }

    const_reference
    operator()() const;

    template <typename... Inds>
    const_reference
    operator()(Inds&&... ii) const;

    reference
    operator()();

    template <typename... Inds>
    reference
    operator()(Inds&&... ii);

    template<typename Indices>
    reference
    operator()(Indices const& ii);

    template<typename Indices>
    const_reference
    operator()(Indices const& ii) const;

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

    storage_type &
    store() { return data_; }

    storage_type const&
    store() const { return data_; }

    private:

    void
    init()
        {
        auto len = area(range_);
#ifdef DEBUG
        if(!isContiguous(range_))
            throw std::runtime_error("Tensor can only be constructed from contiguous range");
        if(len == 0) 
            throw std::runtime_error("Zero area in tensor");
#endif
        data_.assign(len,0.);
        }

    template<typename R>
    void
    assignFromRef(TenRef<const value_type,R> const& ref)
        {
        //TODO: use (optimized) TenRef &= instead ?
        data_.assign(ref.begin(),ref.end());
        }

    void
    assignFrom(Ten const& other)
        {
        if(&other == this) return;
        range_ = other.range_;
        data_ = other.data_;
        }

    void
    moveFrom(Ten && other)
        {
        range_ = std::move(other.range_);
        data_ = std::move(other.data_);
        other.clear();
        }
    };

//
// makeRef functions
//

template<typename T, typename R>
auto
makeRef(TenRef<T,R> & t) { return t; }

//template<typename T, typename R>
//auto
//makeRef(TenRef<const T,R> & t) { return t; }

template<typename T, typename R>
auto
makeRef(Ten<T,R> & t) { return TenRef<T,R>(t); }

template<typename T, typename R>
auto
makeRef(Ten<T,R> const& t) { return TenRef<const T,R>(t); }

//This version of makeRef intended to fail,
//forbids explicitly making TenRefs to temporaries
template<typename T, typename R, typename... Rest>
auto
makeRef(Ten<T,R> && t, Rest&&... args) { return TenRef<T,R>(std::move(t)); }

template<typename T, typename R>
auto
makeRefc(TenRef<T,R> & t) { return TenRef<const T,R>(t); }

template<typename T, typename R>
auto
makeRefc(TenRef<const T,R> & t) { return t; }

template<typename T, typename R>
auto
makeRefc(Ten<T,R> & t) { return TenRef<const T,R>(t); }

template<typename T, typename R>
auto
makeRefc(Ten<T,R> const& t) { return TenRef<const T,R>(t); }

//This version of makeRefc intended to fail,
//forbids explicitly making TenRefs to temporaries
template<typename T, typename R, typename... Rest>
auto
makeRefc(Ten<T,R> && t, Rest&&... args) { return TenRef<const T,R>(std::move(t)); }

template<typename V, typename R>
std::ostream&
operator<<(std::ostream & s, TenRef<V,R> const& T);

template<typename V, typename R>
std::ostream&
operator<<(std::ostream & s, Ten<V,R> const& T);

} //namespace itensor

#include "ten.ih"

#endif
