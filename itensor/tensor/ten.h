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


template<typename range_type>
class Ten;

template<typename range_type>
class TenRefc;

template<typename range_type>
class TenRef;

//Specialize to range_type==Range
using Tensor     = Ten<Range>;
using TensorRef  = TenRef<Range>;
using TensorRefc = TenRefc<Range>;

template<typename Ten_>
using ref_type = typename std::decay<Ten_>::type::ref_type;

template<typename range_type_>
class TenRefc
    { 
    public:
    using range_type = std::remove_const_t<range_type_>;
    using value_type = Real;
    using iterator = TenIter<const Real*,range_type>;
    using const_iterator = iterator;
    using pointer = const Real*;
    using reference = const Real&;
    using size_type = size_t;
    using tensor_type = Ten<range_type>;
    using storage_type = DataRange<const Real>;
    using ref_type = TenRefc<Range>;
    using const_ref_type = TenRefc<Range>;
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

    template<typename Indices>
    reference
    operator()(Indices const& ii) const;

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

template<typename range_type_>
class TenRef : public TenRefc<range_type_>
    { 
    public:
    using parent = TenRefc<range_type_>;
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
    using ref_type = TenRef<Range>;
    using const_ref_type = TenRefc<Range>;

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
    operator=(TenRef const& t) { parent::operator=(t); return *this; }

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
    store() const { return Data(data(),parent::store().size()); }

    pointer
    data() const { return const_cast<pointer>(parent::data()); }

    reference
    operator()() const { return const_cast<reference>(parent::operator()()); }

    template <typename... Inds>
    reference
    operator()(Inds&&... ii) const { return const_cast<reference>(parent::operator()(std::forward<Inds>(ii)...)); }

    template<typename Indices>
    reference
    operator()(Indices const& ii) const { return const_cast<reference>(parent::operator()(ii)); }

    iterator
    begin() const { return iterator(store(),parent::range()); }

    iterator
    end() const { return iterator::makeEnd(parent::range()); }

    };

//Assign to referenced data
template<typename R1, typename R2>
void 
operator&=(TenRef<R1> const& a, TenRefc<R2> const& b);

//Assign to referenced data
template<typename R>
void
operator&=(TenRef<R> const& a, Tensor const& t);

template<typename R1, typename R2>
void 
operator+=(TenRef<R1> const& a, TenRefc<R2> const& b);

template<typename R>
void
operator+=(TenRef<R> const& a, Tensor const& b);

template<typename R>
void
operator+=(Tensor & a, TenRefc<R> const& b);

void
operator+=(Tensor & a, Tensor const& b);

template<typename R1, typename R2, typename Op>
void
transform(TenRefc<R1>  const& from, 
          TenRef<R2> const& to,
          Op&& op);

template<typename range_type>
auto
makeTenRef(Real* p,
           size_t max_size,
           const range_type* prange)
    {
    using R = std::decay_t<std::remove_pointer_t<range_type>>;
    return TenRef<R>({p,max_size},prange);
    }

template<typename range_type>
auto
makeTenRef(const Real* p,
           size_t max_size,
           const range_type* prange)
    {
    using R = std::decay_t<std::remove_pointer_t<range_type>>;
    return TenRefc<R>({p,max_size},prange);
    }

template<typename range_type,
         class = std::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(Real* p,
           size_t max_size,
           range_type && range)
    {
    using R = std::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range type is of pointer type");
    return TenRef<R>({p,max_size},std::move(range));
    }

template<typename range_type,
         class = std::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(const Real* p,
           size_t max_size,
           range_type && range)
    {
    using R = std::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRefc<R>({p,max_size},std::move(range));
    }

template<typename range_type>
auto
makeTenRef(Real* p,
           size_t offset,
           size_t max_size,
           const range_type* prange)
    {
    using R = std::decay_t<std::remove_pointer_t<range_type>>;
    return TenRef<R>({p,offset,max_size},prange);
    }

template<typename range_type>
auto
makeTenRef(const Real* p,
           size_t offset,
           size_t max_size,
           const range_type* prange)
    {
    using R = std::decay_t<std::remove_pointer_t<range_type>>;
    return TenRefc<R>({p,offset,max_size},prange);
    }

template<typename range_type,
         class = std::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(Real* p,
           size_t offset,
           size_t max_size,
           range_type && range)
    {
    using R = std::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRef<R>({p,offset,max_size},std::move(range));
    }

template<typename range_type,
         class = std::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeTenRef(const Real* p,
           size_t offset,
           size_t max_size,
           range_type && range)
    {
    using R = std::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRefc<R>({p,offset,max_size},std::move(range));
    }

template<typename range_type>
auto
makeRef(Data const& store,
        const range_type* prange)
    {
    using R = std::decay_t<std::remove_pointer_t<range_type>>;
    return TenRef<R>(store,prange);
    }

template<typename range_type>
auto
makeRef(cData const& store,
        const range_type* prange)
    {
    using R = std::decay_t<std::remove_pointer_t<range_type>>;
    return TenRefc<R>(store,prange);
    }

template<typename range_type,
         class = std::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeRef(Data const& store,
        range_type && range)
    {
    using R = std::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRef<R>(store,std::move(range));
    }

template<typename range_type,
         class = std::enable_if_t<std::is_rvalue_reference<range_type&&>::value
                               && !std::is_pointer<range_type>::value> >
auto
makeRef(cData const& store,
        range_type && range)
    {
    using R = std::decay_t<range_type>;
    static_assert(!std::is_pointer<R>::value,"Error: range_type is of pointer type");
    return TenRefc<R>(store,std::move(range));
    }

template<typename range_type_>
class Ten
    {
    public:
    using storage_type = std::vector<Real>;
    using ref_storage_type = DataRange<Real>;
    using const_ref_storage_type = DataRange<const Real>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;
    using value_type = Real;
    using pointer = std::add_pointer_t<value_type>;
    using const_pointer = std::add_pointer_t<const value_type>;
    using reference = typename storage_type::reference;
    using const_reference = typename storage_type::const_reference;
    using range_type = range_type_;
    using ref_type = TenRef<range_type>;
    using const_ref_type = TenRefc<range_type>;
    using size_type = typename ref_type::size_type;
    public:
    range_type range_;
    storage_type data_;
    public:

    Ten() { }

    template<typename... Dims>
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
    Ten(TenRefc<R> const& ref) { assignFromRef(ref); }

    template<typename R>
    Ten&
    operator=(TenRefc<R> const& ref) { assignFromRef(ref); return *this; }

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
    assignFromRef(TenRefc<R> const& ref)
        {
        range_ = normalRange(ref.range());
        data_.resize(area(range_));
        makeRef(*this) &= ref;
        }
    };

template<typename R, typename... VArgs>
auto
offset(TenRefc<R> const& t, VArgs&&... vargs) -> decltype(offset(t.range(),0))
    { return offset(t.range(),std::forward<VArgs>(vargs)...); }

template<typename R, typename... VArgs>
auto
offset(Ten<R> const& t, VArgs&&... vargs) -> decltype(offset(t.range(),0))
    { return offset(t.range(),std::forward<VArgs>(vargs)...); }

//
// makeRef functions
//

template<typename R>
auto constexpr
makeRef(TenRef<R> const& t) { return t; }

template<typename R>
auto constexpr
makeRef(TenRefc<R> const& t) { return t; }

template<typename R>
auto 
makeRef(Ten<R> & t) { return TenRef<R>{t}; }

template<typename R>
auto 
makeRef(Ten<R> const& t) { return TenRefc<R>{t}; }

template<typename R, typename Arg, typename... Rest>
auto
makeRef(TenRef<R> const& t, Arg&& arg, Rest&&... rest) 
    { 
    return TenRef<R>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

template<typename R, typename Arg, typename... Rest>
auto
makeRef(TenRefc<R> const& t, Arg&& arg, Rest&&... rest) 
    { 
    return TenRefc<R>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

template<typename R, typename Arg, typename... Rest>
auto
makeRef(Ten<R> & t, Arg&& arg, Rest&&... rest) 
    { 
    return TenRef<R>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

template<typename R, typename Arg, typename... Rest>
auto
makeRef(Ten<R> const& t, Arg&& arg, Rest&&... rest) 
    { 
    return TenRefc<R>(t.store(),std::forward<Arg>(arg),std::forward<Rest>(rest)...); 
    }

//This version of makeRef intended to fail instantiation,
//forbids explicitly making TenRefs to temporaries
template<typename R, typename... VArgs>
auto
makeRef(Ten<R> && t, VArgs&&... args) 
    { 
    static_assert(stdx::false_regardless_of<R,VArgs...>::value,"Cannot call makeRef on temporary/rvalue Ten<R>");
    return TenRefc<R>{};
    }


//
// makeRefc functions
//

template<typename R>
auto
makeRefc(TenRef<R> const& t) { return TenRefc<R>{t}; }

template<typename R>
auto constexpr
makeRefc(TenRefc<R> const& t) { return t; }

template<typename R>
auto
makeRefc(Ten<R> const& t) { return TenRefc<R>{t}; }

//This version of makeRefc intended to fail instantiation,
//forbids explicitly making TenRefs to temporaries
template<typename R, typename... VArgs>
auto
makeRefc(Ten<R> && t, VArgs&&... args) 
    { 
    static_assert(stdx::false_regardless_of<R,VArgs...>::value,"Cannot call makeRefc on temporary/rvalue Ten<R>");
    return TenRefc<R>{};
    }

//
// Other functions
//

template<typename R>
Real
norm(TenRefc<R> const& t);

template<typename R>
Real
norm(Ten<R> const& t) { return norm(makeRefc(t)); }

template<typename R>
bool
isContiguous(TenRefc<R> const& t) { return isContiguous(t.range()); }

template<typename R>
bool
isContiguous(Ten<R> const& t) { return isContiguous(t.range()); }

//Make a scalar (rank 0) tensor with value val
Tensor
scalarTen(Real val);

template<typename R>
void
randomize(TenRef<R> const& t);

template<typename R>
void
randomize(Ten<R> & t);


template<typename R>
std::ostream&
operator<<(std::ostream & s, TenRef<R> const& T);

template<typename R>
std::ostream&
operator<<(std::ostream & s, TenRefc<R> const& T);

inline std::ostream&
operator<<(std::ostream & s, Tensor const& T);

} //namespace itensor

#include "ten.ih"

#endif
