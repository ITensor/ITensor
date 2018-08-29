//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQINDEX_H
#define __ITENSOR_IQINDEX_H
#include "itensor/index.h"
#include "itensor/qn.h"

namespace itensor {

// Forward declarations
struct IndexQN;
class IQIndexDat;
class IQIndexVal;
class IQIndexIter;

template<typename IndexT>
class ITensorT;
using ITensor = ITensorT<Index>;

//
// IQIndex
//

class IQIndex : public Index
    {
    public:
    using storage = std::vector<IndexQN>;
    using storage_ptr = std::shared_ptr<IQIndexDat>;
    using indexval_type = IQIndexVal;
    using const_iterator = IQIndexIter;
    using parent = Index;
    private:
    storage_ptr pd;
    Arrow dir_ = Neither;
    public:

    IQIndex() { }

    //
    // Construct IQIndex from name,
    // collection of any number of 
    // Index-QN pairs.
    // Optional last argument can
    // be Arrow direction of the IQIndex:
    // IQIndex("name",i1,q1,i2,q2,In);
    //
    template<typename... Rest>
    IQIndex(std::string const& name, 
            Index const& i1, QN const& q1, 
            Rest const&... etc);

    // Constructor taking a container
    // of IndexQN's
    IQIndex(std::string const& name, 
            storage && ind_qn, 
            Arrow dir = Out, 
            int plev = 0);


    //number of quantum number blocks
    long 
    nblock() const;

    long 
    nindex() const { return nblock(); }

    //1-indexed
    Index 
    index(long i) const;

    //0-indexed
    Index 
    operator[](long i) const;

    //1-indexed
    QN const& 
    qn(long i) const;

    Arrow 
    dir() const { return dir_; }
    void
    dir(Arrow ndir) { dir_ = ndir; }

    IQIndexVal 
    operator()(long n) const;

    IQIndex& 
    dag();

    const_iterator
    begin() const;

    const_iterator
    end() const;

    void 
    write(std::ostream& s) const;

    IQIndex& 
    read(std::istream& s);

    storage_ptr const&
    store() const { return pd; }

    private:

    void
    makeStorage(storage && iq);

    public:

    //
    // Advanced / developer methods.
    // Not intended for normal usage.
    //

    // Constructor taking a storage pointer
    IQIndex(storage_ptr const& p,
            std::string const& name, 
            Arrow dir = Out, 
            int plev = 0);

    }; //class IQIndex

//
// IndexQN
//

struct IndexQN
    {
    Index index;
    QN qn;

    IndexQN() { }

    IndexQN(Index const& i, 
            QN const& q) 
        : index(i), qn(q) 
        { }

    explicit operator Index() const { return index; }

    void
    dag() { index.dag(); }

    auto
    m() const -> decltype(index.m()) { return index.m(); }

    IndexType
    type() const { return index.type(); }

    void
    write(std::ostream & s) const;

    void
    read(std::istream & s);

    };

bool inline
operator==(IndexQN const& iq, Index const& i) { return iq.index == i; }

bool inline
operator==(Index const& i, IndexQN const& iq) { return iq.index == i; }

bool inline
operator!=(IndexQN const& iq, Index const& i) { return iq.index != i; }

bool inline
operator!=(Index const& i, IndexQN const& iq) { return iq.index != i; }


//
// IQIndexVal
//

class IQIndexVal
    {
    public:
    using index_type = IQIndex;

    IQIndex index;
    long val;

    IQIndexVal();

    IQIndexVal(IQIndex const& iqindex, long val_);

    QN const&
    qn() const;

    explicit operator IndexVal() const;

    explicit operator bool() const { return bool(index); }

    IndexQN
    indexqn() const;

    IndexVal 
    blockIndexVal() const;

    IQIndexVal& 
    dag();

    IQIndexVal& 
    prime(int inc = 1);

    IQIndexVal& 
    prime(IndexType type, int inc = 1);

    IQIndexVal& 
    noprime(IndexType type = All);

    IQIndexVal& 
    mapprime(int plevold, int plevnew, IndexType type = All);
    };

ITensor
operator*(IQIndexVal const& iqiv, IndexVal const& iv);

bool
operator==(IQIndexVal const& iv1, IQIndexVal const& iv2);

bool inline
operator!=(IQIndexVal const& iv1, IQIndexVal const& iv2) { return !(iv1==iv2); }

bool inline
operator==(IQIndexVal const& iv, IQIndex const& I) { return iv.index == I; }

bool inline
operator!=(IQIndexVal const& iv, IQIndex const& I) { return !(iv==I); }

bool inline
operator==(IQIndex const& I, IQIndexVal const& iv) { return iv == I; }

bool inline
operator!=(IQIndex const& I, IQIndexVal const& iv) { return !(iv==I); }

IQIndex inline
dag(IQIndex res) { res.dag(); return res; }

IndexQN inline
dag(IndexQN res) { res.dag(); return res; }

IQIndexVal inline
dag(IQIndexVal res) { res.dag(); return res; }

bool
hasindex(IQIndex const& I, Index const& i);

long
findindex(IQIndex const& I, Index const& i);

long 
offset(IQIndex const& I, Index const& i);

QN 
qn(IQIndex const& I, Index const& i);

Index
findByQN(IQIndex const& I, QN const& qn);

//Make a new index with same properties as I,
//but a different id number (will not compare equal)
//and primelevel zero (or specified value)
//For efficiency, internal sector Index objects
//are the same as I.
IQIndex
sim(IQIndex const& I, int plev = 0);

std::string 
showm(IQIndex const& I);

std::ostream& 
operator<<(std::ostream &o, IQIndex const& I);

std::ostream& 
operator<<(std::ostream &s, IndexQN const& x);

std::ostream& 
operator<<(std::ostream& s, IQIndexVal const& iv);

template<typename... VArgs>
IQIndex
prime(IQIndex I, VArgs&&... vargs) 
    { 
    I.prime(std::forward<VArgs>(vargs)...); 
    return I; 
    }

template<typename... VArgs>
IQIndex
noprime(IQIndex I, VArgs&&... vargs) 
    { 
    I.noprime(std::forward<VArgs>(vargs)...); 
    return I; 
    }

//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
IQIndex inline
mapprime(IQIndex I, 
         int plevold, 
         int plevnew, 
         IndexType type = All)
    { 
    I.mapprime(plevold,plevnew,type); 
    return I; 
    }

template<typename... VArgs>
IQIndexVal
prime(IQIndexVal I, VArgs&&... vargs) 
    { 
    I.prime(std::forward<VArgs>(vargs)...); 
    return I; 
    }

template<typename... VArgs>
IQIndexVal
noprime(IQIndexVal I, VArgs&&... vargs) 
    { 
    I.noprime(std::forward<VArgs>(vargs)...); 
    return I; 
    }

//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
IQIndexVal inline
mapprime(IQIndexVal I, int plevold, int plevnew, IndexType type = All)
    { 
    I.mapprime(plevold,plevnew,type); 
    return I; 
    }

namespace detail {

struct ArrowM
    {
    Arrow dir = Neither;
    long m = 0l;
    ArrowM(Arrow d, long m_) : dir(d), m(m_) { }
    };

ArrowM inline
fill(std::vector<IndexQN> const& v,
     Arrow dir = Out) 
    { 
    return ArrowM(dir,0l);
    }

template<typename... Rest>
ArrowM
fill(std::vector<IndexQN> & v,
     Index const& i, 
     QN const& q, 
     Rest const&... rest)
    {
    v.emplace_back(i,q);
    auto am = fill(v,rest...);
    am.m += i.m();
    return am;
    }

} //namespace detail


template<typename... Rest>
IQIndex::
IQIndex(std::string const& name, 
        Index const& i1, 
        QN const& q1, 
        Rest const&... rest)
    { 
    constexpr auto size = 1+sizeof...(rest)/2;
    auto iq = stdx::reserve_vector<IndexQN>(size);
    auto am = detail::fill(iq,i1,q1,rest...);
    dir_ = am.dir;
    auto I = Index(name,am.m,i1.type(),i1.primeLevel());
    parent::operator=(I);
    makeStorage(std::move(iq));
    }

class IQIndexIter
    { 
    public:
    using iterator_category = std::forward_iterator_tag;
    private:
    IQIndex const& I_;
    long n_ = 0;
    public: 

    IQIndexIter(IQIndex const& I) : I_(I), n_(1l) { }
    
    IndexQN
    operator*() { return IndexQN(I_.index(n_),I_.qn(n_)); }  

    IQIndexIter& 
    operator++() { increment(); return *this; } 

    IQIndexIter 
    operator++(int) { auto prev = *this; increment(); return prev; } 

    bool
    operator==(IQIndexIter const& o) const { return (I_ == o.I_) && (n_ == o.n_); }

    bool
    operator!=(IQIndexIter const& o) const { return !operator==(o); }

    private:

    void
    increment() { ++n_; }

    public:
    //For developer use only; for making end iterator
    IQIndexIter static
    makeEnd(IQIndex const& I)
        {
        auto end = IQIndexIter(I);
        end.n_ = 1+I.nindex();
        return end;
        }
    }; 

IQIndex::const_iterator inline IQIndex::
begin() const { return IQIndexIter(*this); }

IQIndex::const_iterator inline IQIndex::
end() const { return const_iterator::makeEnd(*this); }

} //namespace itensor

#endif
