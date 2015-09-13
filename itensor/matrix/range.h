//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_RANGE_H
#define __ITENSOR_RANGE_H

#include "itensor/util/stdx.h"
#include "itensor/util/readwrite.h"
#include "itensor/tensor/rangeiter.h"

namespace itensor {

template<typename extent_type>
class RangeT;

template<typename extent_type>
class RangeBuilderT;

using Range = RangeT<size_t>;
using RangeBuilder = RangeBuilderT<Range>;

//Storage type for RangeT
template<typename index_type>
struct IndStr
    {
    using size_type = size_t;
    index_type  ind = index_type{}; //convertible to size_type
    size_type   str = 0; //stride

    IndStr(index_type i, size_type s) : ind(i), str(s) { }

    IndStr() { }

    void
    write(std::ostream& s) const
        {
        itensor::write(s,ind);
        itensor::write(s,str);
        }
    void
    read(std::istream& s)
        {
        itensor::read(s,ind);
        itensor::read(s,str);
        }
    };

template<typename index_type_>
class RangeT
    {
    public:
    using index_type = index_type_;
    using value_type = IndStr<index_type>;
    using size_type = typename value_type::size_type;
    using storage_type = InfArray<value_type,8ul>;
    using iterator = RangeIter<RangeT>;
    using const_iterator = iterator;
    private:
    storage_type store_;
    public:

    RangeT() { }

    template <typename... Inds>
    explicit
    RangeT(index_type i0, 
           Inds&&... rest)
        {
        std::array<index_type,1+sizeof...(rest)> inds = 
            {{i0,static_cast<index_type>(rest)...}};
        init(inds);
        }

    explicit
    RangeT(std::vector<index_type> const& v)
        {
        init(v);
        }

    template<size_t size>
    explicit
    RangeT(std::array<index_type,size> const& a)
        {
        init(a);
        }

    explicit
    RangeT(storage_type && store)
      : store_(std::move(store))
        { }

    RangeT(std::initializer_list<index_type> ii) { init(ii); }

    //1-indexed
    size_type
    extent(size_type i) const { return static_cast<size_type>(store_[i-1].ind); }

    //1-indexed
    size_type
    stride(size_type i) const { return store_[i-1].str; }

    size_type
    r() const { return store_.size(); }

    //0-indexed
    value_type const&
    operator[](size_type i) const { return store_[i]; }

    //0-indexed
    value_type &
    operator[](size_type i) { return store_[i]; }

    value_type const&
    at(size_type i) const { return store_.at(i); }

    value_type &
    at(size_type i) { return store_.at(i); }

    size_type
    size() const { return store_.size(); }

    value_type const&
    front() const { return store_.front(); }

    value_type const&
    back() const { return store_.back(); }

    void
    clear() { return store_.clear(); }

    bool
    empty() const { return store_.empty(); }

    value_type*
    data() { return store_.data(); }

    const value_type*
    data() const { return store_.data(); }

    storage_type const&
    store() const { return store_; }

    void
    swap(RangeT& other) { store_.swap(other.store_); }

    template<typename Indexable>
    void 
    init(Indexable const& v);

    void 
    init(std::initializer_list<index_type> il);

    void
    resize(size_type nsize) { store_.resize(nsize); }

    iterator
    begin() const { return RangeIter<RangeT>(*this); }

    iterator
    end() const { return RangeIter<RangeT>::makeEnd(*this); }

    //Compute strides from extents
    void 
    computeStrides()
        {
        size_type str = 1;
        for(auto& i : store_)
            {
            i.str = str;
            str *= static_cast<size_type>(i.ind);
            }
        }
    };

template<typename index_type>
template<typename Indexable>
void RangeT<index_type>::
init(Indexable const& v)
    {
    store_.resize(v.size());
    size_type str = 1;
    for(decltype(v.size()) i = 0; i < v.size(); ++i)
        {
        store_[i].ind = v[i];
        store_[i].str = str;
        str *= static_cast<size_type>(v[i]);
        }
    }

template<typename index_type>
void RangeT<index_type>::
init(std::initializer_list<index_type> il)
    {
    store_.resize(il.size());
    size_type str = 1;
    size_type i = 0;
    for(auto& el : il)
        {
        store_[i].ind = el;
        store_[i].str = str;
        str *= static_cast<size_type>(el);
        ++i;
        }
    }


template<typename index_type>
std::ostream&
operator<<(std::ostream& s, RangeT<index_type> const& r)
    {
    s << "exts: ";
    for(decltype(r.r()) i = 1; i <= r.r(); ++i) s << r.extent(i) << " ";
    s << "strs: ";
    for(decltype(r.r()) i = 1; i <= r.r(); ++i) s << r.stride(i) << " ";
    return s;
    }


template<typename index_type>
void
write(std::ostream& s, RangeT<index_type> const& r)
    {
    itensor::write(s,r.store());
    }

template<typename index_type>
void
read(std::istream& s, RangeT<index_type>& r)
    {
    using storage_type = typename RangeT<index_type>::storage_type;
    storage_type store;
    itensor::read(s,store);
    r = RangeT<index_type>(std::move(store));
    }


template<typename range_type_>
class RangeBuilderT
    {
    public:
    using range_type = typename std::decay<range_type_>::type;
    using index_type = typename range_type::index_type;
    using size_type = typename range_type::size_type;
    using storage_type = typename range_type::storage_type;
    using value_type = typename storage_type::value_type;
    private:
    storage_type store_;
    size_t n_ = 0;
    bool auto_compute_strides_ = true;
    public:

    RangeBuilderT() { }

    RangeBuilderT(size_type size)
      : store_(size)
        { }

    //"build" and return the range
    //store_ will be moved into the newly created range
    range_type
    build()
        {
        auto res = range_type{std::move(store_)};
        if(auto_compute_strides_) res.computeStrides();
        return res;
        }

    operator bool() const { return !store_.empty(); }

    size_type
    size() const { return store_.size(); }

    void
    resize(size_type newsize) { store_.resize(newsize); }

    //1-indexed
    size_type const&
    extent(size_type j) const { return static_cast<size_type>(store_[j-1].ind); }

    //1-indexed
    size_type
    stride(size_type j) const { return store_[j-1].str; }

    storage_type &
    store() { return store_; }

    storage_type const&
    store() const { return store_; }

    void
    nextIndex(index_type const& ind)
        {
        store_[n_].ind = ind;
        ++n_;
        }

    void
    nextIndStr(index_type const& ind,
               size_type str)
        {
        store_[n_].ind = ind;
        store_[n_].str = str;
        auto_compute_strides_ = false;
        ++n_;
        }

    void
    setIndex(size_type j, index_type const& i)
        {
        store_[j].ind = i;
        }

    void
    setStride(size_type j, size_type s)
        {
        store_[j].str = s;
        auto_compute_strides_ = false;
        }

    void
    setIndStr(size_type j, 
              index_type const& ind,
              size_type str)
        {
        store_[j].ind = ind;
        store_[j].str = str;
        auto_compute_strides_ = false;
        }

    void
    sortByIndex()
        {
        auto comp = [](value_type const& i1, value_type const& i2) { return i1.ind > i2.ind; };
        std::sort(store_.begin(),store_.end(),comp);
        auto_compute_strides_ = true;
        }
    };

namespace detail {

    template<typename RangeT, typename Iterable>
    auto
    offsetIterable(RangeT const& r, Iterable const& inds)
        {
        using size_type = typename RangeT::size_type;
        size_type I  = 0, 
                  ri = 1;
        for(auto& ii : inds)
            {
#ifdef DEBUG
            if(ri > size_type(r.r()))
                Error("Container-Range size mismatch in offset(...)");
#endif
            I += r.stride(ri) * (ii-1);
            ++ri;
            }
        return I;
        }

    template<typename RangeT>
    struct ComputeOffset
        {
        using size_type = typename RangeT::size_type;

        RangeT const& r;
        ComputeOffset(RangeT const& r_) : r(r_) { }

        template<typename... Inds>
        size_type
        operator()(size_type first, Inds... rest) const 
            { 
            return off<0>(first,rest...);
            }

        private:

        template <size_type i, typename... Inds>
        size_type
        off(size_type first, Inds... rest) const
            {
            return (first-1)*r.stride(i) + off<i+1>(rest...);
            }

        template <size_type i>
        size_type
        off(size_type first) const
            {
            return (first-1)*r.stride(i);
            }
        };

} //namespace detail


//1-indexed
template<typename RangeT, typename Iterable>
auto
offset(RangeT const& r, Iterable const& inds)
    //Constrain this template to only work for inds that have a begin() method
    -> stdx::if_compiles_return<decltype(inds.begin()),decltype(r.extent(1))>
    //...if so make the return type to be decltype(r.extent(1))
    {
    return detail::offsetIterable(r,inds);
    }

//1-indexed
template<typename RangeT, typename... Inds>
auto
offset(RangeT const& r, size_t i1, Inds... inds)
    {
    return detail::ComputeOffset<RangeT>(r)(i1,inds...);
    }

template<typename RangeT>
auto
area(RangeT const& R)
    { 
    using size_type = typename RangeT::size_type;
    size_type A = 1;
    for(decltype(R.r()) n = 1; n <= R.r(); ++n)
        {
        A *= R.extent(n);
        }
    return A;
    }

//make Range with same extents but
//normal (unsliced) strides
template<typename index_type>
Range
normalRange(RangeT<index_type> const& R)
    {
    auto rb = RangeBuilder(R.r());
    for(decltype(R.r()) n = 1; n <= R.r(); ++n)
        rb.nextIndex(R.extent(n));
    return rb.build();
    }


//
//A range R is contiguous if collecting
//all possible outputs of offset(R,...) yields
//the set {0,1,...,area(R)-1} (though 
//in no particular order)
//For this to be true, sufficient that
//the max possible output of offset(R,...)
//equals (area(R)-1)
//
//Proof:
// Ranges always start at offset(R,{1,1,1...})=0;
// area(R) gives the number of outputs;
// IF the max output is area(R)-1 and
// there are area(R) outputs, the only
// set fulfilling this is {0,1,...,area(R)-1})
//
template<typename RangeT>
bool
isContiguous(RangeT const& R)
    {
    using size_type = typename RangeT::size_type;
    size_type max_offset = 0,
              area = 1;
    for(decltype(R.r()) n = 1; n <= R.r(); ++n)
        {
        max_offset += R.stride(n)*(R.extent(n)-1);
        area *= R.extent(n);
        }
    return (1+max_offset) == area;
    }


} //namespace itensor

#endif
