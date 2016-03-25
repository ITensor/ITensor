//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQINDEX_H
#define __ITENSOR_IQINDEX_H
#include "itensor/index.h"
#include "itensor/qn.h"
#include "itensor/util/readwrite.h"

namespace itensor {

// Forward declarations
struct IndexQN;
class IQIndexDat;
class IQIndexVal;
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
    using const_iterator = storage::const_iterator;
    private:
    storage_ptr pd;
    Arrow dir_;
    public:

    IQIndex();

    //
    //Construct IQIndex from name,
    //collection of any number of 
    //Index-QN pairs.
    //Optional last argument can
    //be Arrow direction of the IQIndex:
    //IQIndex("name",i1,q1,i2,q2,In);
    //
    template<typename... Rest>
    IQIndex(std::string const& name, 
            Index const& i1, QN const& q1, 
            Rest const&... etc);

    //Constructor taking a container
    //of IndexQN's
    IQIndex(std::string const& name, 
            storage && ind_qn, 
            Arrow dir = Out, 
            int plev = 0);

    //
    //Accessor Methods
    //

    long 
    nindex() const;

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

    //
    // Operators
    //

    IQIndexVal 
    operator()(long n) const;

    //
    // Other methods
    //

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

    }; //class IQIndex

//
// IndexQN
//

struct IndexQN //: public Index
    {
    using IndexValT = IQIndexVal;

    Index index;
    QN qn;

    IndexQN() { }

    IndexQN(Index const& i, QN const& q) : index(i), qn(q) { }

    explicit operator Index() const { return index; }

    void
    dag() { index.dag(); }

    auto
    m() const -> decltype(index.m()) { return index.m(); }

    IndexType
    type() const { return index.type(); }

    void 
    write(std::ostream& s) const { index.write(s); itensor::write(s,qn); }

    void 
    read(std::istream& s) { index.read(s); itensor::read(s,qn); }
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

    explicit operator IndexVal() const;

    explicit operator bool() const { return bool(index); }

    QN const&
    qn() const;

    QN const&
    qn(long j) const { return index.qn(j); }

    IndexQN
    indexqn() const;

    IndexVal 
    blockIndexVal() const;

    long 
    m() const { return index.m(); }

    IQIndexVal& 
    prime(int inc = 1);

    IQIndexVal& 
    prime(IndexType type, int inc = 1);

    IQIndexVal& 
    noprime(IndexType type = All);

    IQIndexVal& 
    mapprime(int plevold, int plevnew, IndexType type = All);

    IQIndexVal& 
    dag();

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
hasindex(const IQIndex& I, const Index& i);

long
findindex(const IQIndex& I, const Index& i);

long 
offset(const IQIndex& I, const Index& i);

QN 
qn(const IQIndex& I, const Index& i);

Index
findByQN(const IQIndex& I, const QN& qn);

std::string 
showm(IQIndex const& I);

std::ostream& 
operator<<(std::ostream &o, const IQIndex &I);

std::ostream& 
operator<<(std::ostream &s, const IndexQN& x);

std::ostream& 
operator<<(std::ostream& s, const IQIndexVal& iv);

template<typename... VarArgs>
IQIndex
prime(IQIndex I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
IQIndex
noprime(IQIndex I, VarArgs&&... vargs) { I.noprime(std::forward<VarArgs>(vargs)...); return I; }

//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
IQIndex inline
mapprime(IQIndex I, int plevold, int plevnew, IndexType type = All)
    { I.mapprime(plevold,plevnew,type); return I; }

template<typename... VarArgs>
IQIndexVal
prime(IQIndexVal I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
IQIndexVal
noprime(IQIndexVal I, VarArgs&&... vargs) { I.noprime(std::forward<VarArgs>(vargs)...); return I; }

//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
IQIndexVal inline
mapprime(IQIndexVal I, int plevold, int plevnew, IndexType type = All)
    { I.mapprime(plevold,plevnew,type); return I; }


//
//
// Implementations
//
//

class IQIndexDat
    {
    public:
    using storage = std::vector<IndexQN>;
    using iterator = storage::iterator;
    using const_iterator = storage::const_iterator;
    private:
    storage iq_;
    public:

    IQIndexDat() { }

    template<typename... Rest>
    IQIndexDat(Index const& i1, QN const& q1,
               Rest const&... args) 
        { 
        constexpr auto size = sizeof...(args)/2+1;
        iq_.resize(size);
        fill<0>(i1,q1,args...);
        }

    explicit
    IQIndexDat(const storage& ind_qn) 
      : iq_(ind_qn)
        { }

    explicit
    IQIndexDat(storage&& ind_qn) 
      : iq_(std::move(ind_qn))
        { }

    //Disallow copying
    IQIndexDat(const IQIndexDat&) = delete;

    void 
    operator=(IQIndexDat const&) = delete;

    storage const&
    inds() const { return iq_; }

    long
    size() { return iq_.size(); }

    Index const&
    index(long i) { return iq_[i-1].index; }

    Index const&
    operator[](long i) { return iq_[i].index; }

    const QN&
    qn(long i) { return iq_[i-1].qn; }

    iterator
    begin() { return iq_.begin(); }

    iterator
    end() { return iq_.end(); }

    const_iterator
    begin() const { return iq_.begin(); }

    const_iterator
    end()   const { return iq_.end(); }

    void 
    write(std::ostream& s) const { itensor::write(s,iq_); }

    void 
    read(std::istream& s) { itensor::read(s,iq_); }

    //IQIndex::storage_ptr
    //clone() { return std::make_shared<IQIndexDat>(iq_); }

    template<long J>
    void
    fill(Arrow dir = Out) { }

    template<long J, typename... Rest>
    void
    fill(Index const& i, 
         QN const& q, 
         Rest const&... rest)
        {
#ifdef DEBUG
        assert(J < iq_.size());
#endif
        iq_[J] = IndexQN(i,q);
        fill<J+1>(rest...);
        }
    };

long
totalM(IQIndexDat::storage const& storage);

inline IQIndex::
IQIndex() 
    : 
    dir_(Neither)
    { }

namespace detail {

long inline
totalM(Arrow dir = Out)
    {
    return 0;
    }

template<typename... Rest>
long
totalM(Index const& i, 
       QN const& q,
       Rest const&... rest)
    {
    return i.m()+detail::totalM(rest...);
    }

Arrow inline
getDir(Arrow dir = Out)
    {
    return dir;
    }

template<typename... Rest>
Arrow
getDir(Index const& i, 
       QN const& q,
       Rest const&... rest)
    {
    return detail::getDir(rest...);
    }

} //namespace detail


template<typename... Rest>
IQIndex::
IQIndex(std::string const& name, 
        Index const& i1, QN const& q1, 
        Rest const&... rest)
  : Index(name,detail::totalM(i1,q1,rest...),i1.type(),i1.primeLevel()), 
    pd(std::make_shared<IQIndexDat>(i1,q1,rest...)),
    dir_(detail::getDir(rest...))
    { }

inline IQIndex::
IQIndex(std::string const& name, 
        storage && ind_qn, 
        Arrow dir, 
        int plev) 
  : Index(name,totalM(ind_qn),ind_qn.front().index.type(),plev),
    pd(std::make_shared<IQIndexDat>(std::move(ind_qn))),
    dir_(dir)
    { }


} //namespace itensor

#endif
