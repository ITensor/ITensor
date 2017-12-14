//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEN_H_
#define __ITENSOR_TEN_H_

#include "itensor/detail/algs.h"
#include "itensor/tensor/teniter.h"
#include "itensor/tensor/range.h"
#include "itensor/tensor/lapack_wrap.h"

namespace itensor {

template<typename range_type, typename T = Real>
class Ten;

template<typename range_type, typename T = Real>
class TenRefc;

template<typename range_type, typename T = Real>
class TenRef;

//Specialize to range_type==Range
using Tensor     = Ten<Range,Real>;
using TensorRef  = TenRef<Range,Real>;
using TensorRefc = TenRefc<Range,Real>;

using CTensor     = Ten<Range,Cplx>;
using CTensorRef  = TenRef<Range,Cplx>;
using CTensorRefc = TenRefc<Range,Cplx>;

//1-indexed versions
using Tensor1     = Ten<Range1,Real>;
using TensorRef1  = TenRef<Range1,Real>;
using TensorRefc1 = TenRefc<Range1,Real>;

using CTensor1     = Ten<Range1,Cplx>;
using CTensorRef1  = TenRef<Range1,Cplx>;
using CTensorRefc1 = TenRefc<Range1,Cplx>;

class TensorType { };

template<typename Derived>
struct isTensor
    {
    constexpr isTensor() { }
    bool static constexpr value = std::is_base_of<TensorType,stdx::decay_t<Derived>>::value;
    constexpr operator bool() const noexcept { return value; }
    };

template<typename Ten_, typename range_type = Range>
using ref_type = typename stdx::decay_t<Ten_>::template ref_type<range_type>;

template<typename Ten_>
using range_type = typename stdx::decay_t<Ten_>::range_type;



template<typename T, bool istensor = isTensor<T>::value >
struct ValTypeHelper { using type = typename stdx::decay_t<T>::value_type; };
template<typename T>
struct ValTypeHelper<T,false> { using type = T; };
template<typename T>
using val_type = typename ValTypeHelper<T>::type;

template<typename TA, typename TB>
using common_type = stdx::conditional_t<(std::is_same<val_type<TA>,Cplx>::value || std::is_same<val_type<TB>,Cplx>::value),
                                        Cplx,
                                        Real>;



template<typename range_type_,typename value_type_>
class TenRefc : public TensorType
    { 
    static_assert(not std::is_const<value_type_>::value,
                  "value type template argument of TenRefc should not be const");
    public:
    using range_type = stdx::remove_const_t<range_type_>;
    using value_type = value_type_;
    using iterator = TenIter<const value_type*,range_type>;
    using const_iterator = iterator;
    using pointer = const value_type*;
    using reference = const value_type&;
    using size_type = size_t;
    using tensor_type = Ten<range_type,value_type>;
    using storage_type = DataRange<const value_type>;
    template<typename R_>
    using ref_type = TenRefc<R_,value_type>;
    template<typename R_>
    using const_ref_type = TenRefc<R_,value_type>;
    private:
    storage_type d_;
    const range_type* prange_ = nullptr;
    range_type range_;
    public:

    TenRefc() { }

    TenRefc(storage_type dat,
            range_type const& range)
      : d_(dat),
        prange_(&range)
        { }

    TenRefc(storage_type dat,
            range_type && range)
      : d_(dat),
        range_(std::move(range))
        { 
        prange_ = &range_;
        }

    TenRefc(storage_type dat,
            const range_type* prange)
      : d_(dat),
        prange_(prange)
        { }

    TenRefc(TenRefc const& t)
        {
        operator=(t);
        }

    TenRefc&
    operator=(TenRefc const& t);

    TenRefc(TenRefc && t)
        {
        operator=(std::move(t));
        }

    TenRefc&
    operator=(TenRefc && t);

    TenRefc(tensor_type const& t) { pointTo(t); }

    TenRefc&
    operator=(tensor_type const& t) { pointTo(t); return *this; }

    TenRefc(tensor_type&& t) = delete;

    TenRefc&
    operator=(tensor_type&& t) = delete;

    bool
    ownRange() const { return prange_ == &range_; }

    size_type
    r() const { return prange_->r(); }

    size_type 
    size() const { return area(*prange_); }

    explicit operator bool() const { return bool(d_.data());}

    size_type
    extent(size_type i) const { return prange_->extent(i); }

    size_type
    stride(size_type i) const { return prange_->stride(i); }

    range_type const&
    range() const { return *prange_; }

    storage_type const&
    store() const { return d_; }

    // direct access to data
    pointer
    data() const { return d_.data(); }

    reference
    operator()() const;

    template <typename... Inds>
    reference
    operator()(Inds&&... ii) const;

    reference
    operator[](size_t n) const;

    void
    clear() { d_.clear(); prange_ = nullptr; }

    const_iterator
    begin() const { return iterator(d_,*prange_); }

    const_iterator
    end() const { return iterator::makeEnd(*prange_); }

    const_iterator
    cbegin() const { return const_iterator(d_,*prange_); }

    const_iterator
    cend() const { return const_iterator::makeEnd(*prange_); }

    void
    pointTo(tensor_type const& t);
    };

template<typename range_type_, typename value_type_>
class TenRef : public TenRefc<range_type_,value_type_>
    { 
    static_assert(not std::is_const<value_type_>::value,
                  "value type template argument of TenRef should not be const");
    public:
    using parent = TenRefc<range_type_,value_type_>;
    using range_type = typename parent::range_type;
    using value_type = typename parent::value_type;
    using iterator = TenIter<value_type*,range_type>;
    using const_iterator = TenIter<const value_type*,range_type>;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = typename parent::size_type;
    using tensor_type = typename parent::tensor_type;
    using storage_type = DataRange<value_type>;
    template<typename R_>
    using ref_type = TenRef<R_,value_type>;
    template<typename R_>
    using const_ref_type = TenRefc<R_,value_type>;

    TenRef() { }

    TenRef(storage_type dat,
           range_type const& range)
      : parent(dat,range)
        { }

    TenRef(storage_type dat,
           range_type && range)
      : parent(dat,std::move(range))
        { }

    TenRef(storage_type dat,
            const range_type* prange)
      : parent(dat,prange)
        { }

    TenRef(TenRef const& t)
        {
        operator=(t);
        }

    TenRef&
    operator=(TenRef const& t) 
        { 
        parent::operator=(t); 
        return *this; 
        }

    TenRef(TenRef && t)
        {
        operator=(std::move(t));
        }

    TenRef&
    operator=(TenRef && t) { parent::operator=(std::move(t)); return *this; }

    TenRef(tensor_type & t) { parent::pointTo(t); }

    TenRef&
    operator=(tensor_type & t) { parent::pointTo(t); return *this; }

    storage_type
    store() const { return storage_type(data(),parent::store().size()); }

    pointer
    data() const { return const_cast<pointer>(parent::data()); }

    reference
    operator()() const { return const_cast<reference>(parent::operator()()); }

    template <typename... Inds>
    reference
    operator()(Inds&&... ii) const 
        { 
        return const_cast<reference>(parent::operator()(std::forward<Inds>(ii)...)); 
        }

    reference
    operator[](size_t n) const;

    iterator
    begin() const { return iterator(store(),parent::range()); }

    iterator
    end() const { return iterator::makeEnd(parent::range()); }

    };

//Assign to referenced data
template<typename R1, typename R2, typename T>
void 
operator&=(TenRef<R1,T> const& A, TenRefc<R2,T> const& B);

//Assign to referenced data
template<typename R1, typename R2, typename T>
void
operator&=(TenRef<R1,T> const& A, Ten<R2,T> const& B);

template<typename R1, typename R2, typename T>
void 
operator+=(TenRef<R1,T> const& a, TenRefc<R2,T> const& b);

template<typename R, typename T>
void
operator+=(TenRef<R,T> const& a, Ten<Range,T> const& b);

template<typename R,typename T>
void
operator+=(Ten<Range,T> & a, TenRefc<R,T> const& b);

//template<typename R, typename T>
//void
//operator+=(Ten<R,T> & a, Ten<R,T> const& b);

template<typename R1, typename T1, 
         typename R2, typename T2, 
         typename Op>
void
transform(TenRefc<R1,T1> const& from, 
          TenRef<R2,T2>  const& to,
          Op&& op);

template<typename V, typename range_type>
auto
makeTenRef(V * p,
           size_t max_size,
           const range_type* prange)
    -> TenRef<stdx::decay_t<stdx::remove_pointer_t<range_type>>,V>
    {
    using R = stdx::decay_t<stdx::remove_pointer_t<range_type>>;
    return TenRef<R,V>({p,max_size},prange);
    }

template<typename V, typename range_type>
auto
makeTenRef(V const* p,
           size_t max_size,
           range_type const* prange)
    -> TenRefc<stdx::decay_t<stdx::remove_pointer_t<range_type>>,V>
    {
    using R = stdx::decay_t<stdx::remove_pointer_t<range_type>>;
    return TenRefc<R,V>({p,max_size},prange);
    }

template<typename V, typename range_type,
         class = stdx::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(V * p,
           size_t max_size,
           range_type && range)
    -> TenRef<stdx::decay_t<range_type>,V>
    {
    using R = stdx::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range type is of pointer type");
    return TenRef<R,V>({p,max_size},std::move(range));
    }

template<typename V, typename range_type,
         class = stdx::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(V const* p,
           size_t max_size,
           range_type && range)
    -> TenRefc<stdx::decay_t<range_type>,V>
    {
    using R = stdx::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRefc<R,V>({p,max_size},std::move(range));
    }

template<typename V, typename range_type>
auto
makeTenRef(V * p,
           size_t offset,
           size_t max_size,
           range_type const* prange)
    -> TenRef<stdx::decay_t<stdx::remove_pointer_t<range_type>>,V>
    {
    using R = stdx::decay_t<stdx::remove_pointer_t<range_type>>;
    return TenRef<R,V>({p,offset,max_size},prange);
    }

template<typename V, typename range_type>
auto
makeTenRef(V const* p,
           size_t offset,
           size_t max_size,
           range_type const* prange)
    -> TenRefc<stdx::decay_t<stdx::remove_pointer_t<range_type>>,V>
    {
    using R = stdx::decay_t<stdx::remove_pointer_t<range_type>>;
    return TenRefc<R,V>({p,offset,max_size},prange);
    }

template<typename V, typename range_type,
         class = stdx::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(V * p,
           size_t offset,
           size_t max_size,
           range_type && range)
    -> TenRef<stdx::decay_t<range_type>,V>
    {
    using R = stdx::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRef<R,V>({p,offset,max_size},std::move(range));
    }

template<typename V, typename range_type,
         class = stdx::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(V const* p,
           size_t offset,
           size_t max_size,
           range_type && range)
    -> TenRefc<stdx::decay_t<range_type>,V>
    {
    using R = stdx::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRefc<R,V>({p,offset,max_size},std::move(range));
    }

template<typename range_type,typename T>
auto
makeRef(DataRange<T> const& store,
        range_type const* prange)
    -> TenRef<stdx::decay_t<stdx::remove_pointer_t<range_type>>,T>
    {
    using R = stdx::decay_t<stdx::remove_pointer_t<range_type>>;
    return TenRef<R,T>(store,prange);
    }

template<typename range_type,typename T>
auto
makeRef(DataRange<const T> const& store,
        range_type const* prange)
    -> TenRefc<stdx::decay_t<stdx::remove_pointer_t<range_type>>,T>
    {
    using R = stdx::decay_t<stdx::remove_pointer_t<range_type>>;
    return TenRefc<R,T>(store,prange);
    }

template<typename range_type,
         typename T,
         class = stdx::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeRef(DataRange<T> const& store,
        range_type && range)
    -> TenRef<stdx::decay_t<range_type>,T>
    {
    using R = stdx::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRef<R,T>(store,std::move(range));
    }

template<typename range_type,
         typename T,
         class = stdx::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeRef(DataRange<const T> const& store,
        range_type && range)
    -> TenRefc<stdx::decay_t<range_type>,T>
    {
    using R = stdx::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRefc<R,T>(store,std::move(range));
    }

template<typename range_type_, typename value_type_>
class Ten : public TensorType
    {
    public:
    using value_type = value_type_;
    using storage_type = std::vector<value_type>;
    using ref_storage_type = DataRange<value_type>;
    using const_ref_storage_type = DataRange<const value_type>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;
    using pointer = stdx::add_pointer_t<value_type>;
    using const_pointer = stdx::add_pointer_t<const value_type>;
    using reference = typename storage_type::reference;
    using const_reference = typename storage_type::const_reference;
    using range_type = range_type_;
    template<typename R_>
    using ref_type = TenRef<R_,value_type>;
    template<typename R_>
    using const_ref_type = TenRefc<R_,value_type>;
    using size_type = typename ref_type<range_type>::size_type;
    private:
    static_assert(not std::is_const<value_type>::value,
                  "value type template argument of Ten should not be const");
    range_type range_;
    storage_type data_;
    public:

    Ten() { }

    template<typename... Dims>
    explicit
    Ten(size_type d1, Dims&&... rest)
      : range_(d1,std::forward<Dims>(rest)...)
        {
        init();
        }

    Ten(storage_type && store,
          range_type && range) 
      : range_(std::move(range)),
        data_(std::move(store))
        { 
#ifdef DEBUG
        if(!isContiguous(range_)) Error("Tensor required to have contiguous range");
#endif
        }

    template<typename R>
    explicit
    Ten(TenRefc<R,value_type> const& ref) { assignFromRef(ref); }

    template<typename R>
    Ten&
    operator=(TenRefc<R,value_type> const& ref) { assignFromRef(ref); return *this; }

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

    value_type
    operator()() const;

    template <typename... Inds>
    const_reference
    operator()(Inds&&... ii) const;

    reference
    operator()();

    template <typename... Inds>
    reference
    operator()(Inds&&... ii);

    const_reference
    operator[](size_t n) const;

    reference
    operator[](size_t n);


    //template<typename Indices>
    //reference
    //operator()(Indices const& ii);

    //template<typename Indices>
    //const_reference
    //operator()(Indices const& ii) const;

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

    ref_storage_type
    store() { return ref_storage_type(data(),size()); }

    const_ref_storage_type
    store() const { return const_ref_storage_type(data(),size()); }

    storage_type &
    storage() { return data_; }

    storage_type const&
    storage() const { return data_; }

    void
    resize(range_type const& newrange)
        {
        range_ = newrange;
        data_.resize(area(range_));
        }

    void
    swap(Ten & other)
        {
        range_.swap(other.range_);
        data_.swap(other.data_);
        }

    void
    read(std::istream& s)
        {
        itensor::read(s,range_);
        itensor::read(s,data_);
        }

    void
    write(std::ostream& s) const
        {
        itensor::write(s,range_);
        itensor::write(s,data_);
        }

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
    assignFromRef(TenRefc<R,value_type> const& ref)
        {
        range_ = normalRange(ref.range());
        data_.resize(area(range_));
        makeRef(*this) &= ref;
        }
    };

template<typename R, typename T, typename... VArgs>
auto
offset(TenRefc<R,T> const& t, VArgs&&... vargs) -> size_t
    { return offset(t.range(),std::forward<VArgs>(vargs)...); }

template<typename R, typename T, typename... VArgs>
auto
offset(Ten<R,T> const& t, VArgs&&... vargs) -> size_t
    { return offset(t.range(),std::forward<VArgs>(vargs)...); }

//
// makeRef functions
//

template<typename R,typename T>
auto constexpr
makeRef(TenRef<R,T> const& t) -> decltype(t) { return t; }

template<typename R,typename T>
auto constexpr
makeRef(TenRefc<R,T> const& t) -> decltype(t) { return t; }

template<typename R,typename T>
auto 
makeRef(Ten<R,T> & t) -> TenRef<R,T> { return TenRef<R,T>{t}; }

template<typename R,typename T>
auto 
makeRef(Ten<R,T> const& t) -> TenRefc<R,T> { return TenRefc<R,T>{t}; }

template<typename R, typename T, typename Arg, typename... Rest>
auto
makeRef(TenRef<R,T> const& t, Arg&& arg, Rest&&... rest) 
    -> TenRef<R,T>
    { 
    return TenRef<R,T>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

template<typename R, typename T, typename Arg, typename... Rest>
auto
makeRef(TenRefc<R,T> const& t, Arg&& arg, Rest&&... rest) 
    -> TenRefc<R,T>
    { 
    return TenRefc<R,T>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

template<typename R,typename T, typename Arg, typename... Rest>
auto
makeRef(Ten<R,T> & t, Arg&& arg, Rest&&... rest) 
    -> TenRef<R,T>
    { 
    return TenRef<R,T>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

template<typename R,typename T, typename Arg, typename... Rest>
auto
makeRef(Ten<R,T> const& t, Arg&& arg, Rest&&... rest) 
    -> TenRefc<R,T>
    { 
    return TenRefc<R,T>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

//This version of makeRef intended to fail instantiation,
//forbids explicitly making TenRefs to temporaries
template<typename R,typename T, typename... VArgs>
auto
makeRef(Ten<R,T> && t, VArgs&&... args) 
    -> TenRefc<R,T>
    { 
    static_assert(stdx::false_regardless_of<R,VArgs...>::value,"Cannot call makeRef on temporary/rvalue Ten object");
    return TenRefc<R,T>{};
    }


//
// makeRefc functions
//

template<typename R,typename T>
auto
makeRefc(TenRef<R,T> const& t) -> TenRefc<R,T> { return TenRefc<R,T>{t}; }

template<typename R, typename T>
auto constexpr
makeRefc(TenRefc<R,T> const& t) -> decltype(t) { return t; }

template<typename R,typename T>
auto
makeRefc(Ten<R,T> const& t) -> TenRefc<R,T> { return TenRefc<R,T>{t}; }

//This version of makeRefc intended to fail instantiation,
//forbids explicitly making TenRefs to temporaries
template<typename R,typename T, typename... VArgs>
auto
makeRefc(Ten<R,T> && t, VArgs&&... args) 
    -> TenRefc<R,T>
    { 
    static_assert(stdx::false_regardless_of<R,T,VArgs...>::value,"Cannot call makeRefc on temporary/rvalue Ten<R>");
    return TenRefc<R,T>{};
    }

//
// Other functions
//

template<typename R,typename T>
auto
rank(TenRefc<R,T> const& t) -> decltype(rank(t.range())) { return rank(t.range()); }

template<typename R,typename T>
auto
rank(Ten<R,T> const& t) -> decltype(rank(t.range())) { return rank(t.range()); }

//ord is alias for rank, order is preferred in applied math literature over rank
template<typename R,typename T>
auto
ord(TenRefc<R,T> const& t) -> decltype(rank(t.range())) { return rank(t.range()); }

template<typename R,typename T>
auto
ord(Ten<R,T> const& t) -> decltype(rank(t.range())) { return rank(t.range()); }

template<typename R, typename V>
Real
norm(TenRefc<R,V> const& t);

template<typename R, typename T>
Real
norm(Ten<R,T> const& t) { return norm(makeRefc(t)); }

template<typename R, typename T>
bool
isContiguous(TenRefc<R,T> const& t) { return isContiguous(t.range()); }

template<typename R, typename T>
bool
isContiguous(Ten<R,T> const& t) { return isContiguous(t.range()); }

//Make a scalar (rank 0) tensor with value val
Tensor
scalarTen(Real val);

template<typename R, typename V>
void
randomize(TenRef<R,V> const& t);

template<typename R, typename V>
void
randomize(Ten<R,V> & t);

template<typename R>
Data
realData(TenRef<R,Real> const& t) { return Data(t.data(),t.size()); }

template<typename R>
Datac
realData(TenRefc<R,Real> const& t) { return Datac(t.data(),t.size()); }

template<typename R>
Data
realData(TenRef<R,Cplx> const& t) { return Data(reinterpret_cast<Real*>(t.data()),2*t.size()); }

template<typename R>
Datac
realData(TenRefc<R,Cplx> const& t) { return Datac(reinterpret_cast<const Real*>(t.data()),2*t.size()); }

template<typename T, class = stdx::require<isTensor<T>> >
bool constexpr
isReal(T const& t) { return std::is_same<typename T::value_type,Real>::value; }

template<typename T, class = stdx::require<isTensor<T>> >
bool constexpr
isCplx(T const& t) { return std::is_same<typename T::value_type,Cplx>::value; }

//conjugate in-place, modifying elements
template<typename R, typename V>
void
conjugate(TenRef<R,V> const& T);

template<typename R, typename V>
void
conjugate(Ten<R,V> & T);

//return conjugated copy
template<typename R, typename V>
Ten<R,V>
conj(TenRefc<R,V> const& T);

template<typename R, typename V>
Ten<R,V>
conj(Ten<R,V> T);

template<typename R, typename V>
std::ostream&
operator<<(std::ostream & s, TenRef<R,V> const& T);

template<typename R, typename V>
std::ostream&
operator<<(std::ostream & s, TenRefc<R,V> const& T);

template<typename R, typename V>
std::ostream&
operator<<(std::ostream & s, Ten<R,V> const& T);

} //namespace itensor

#include "ten.ih"

#endif
