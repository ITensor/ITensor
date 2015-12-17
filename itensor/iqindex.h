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

    template<typename... Rest>
    IQIndex(std::string const& name, 
            Index const& i1, QN const& q1, 
            Rest const&... etc);

    template<typename... Rest>
    IQIndex(std::string const& name, 
            Arrow dir,
            Index const& i1, QN const& q1, 
            Rest const&... etc);

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

    //////////
    IQIndex index;
    long val;
    //////////

    IQIndexVal();

    IQIndexVal(const IQIndex& iqindex, long val_);

    const QN&
    qn() const;

    const QN&
    qn(long j) const { return index.qn(j); }

    bool
    operator==(const IQIndexVal& other) const;
    bool
    operator!=(const IQIndexVal& other) const { return !operator==(other); }

    bool
    operator==(const IQIndex& iqind) const { return index == iqind; }
    bool
    operator!=(const IQIndex& iqind) const { return !operator==(iqind); }

    explicit operator bool() const { return bool(index); }

    IndexQN
    indexqn() const;

    operator IndexVal() const;

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

    ITensor 
    operator*(const IndexVal& iv) const;

    };


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
showm(const IQIndex& I);

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

    template<typename... Args>
    IQIndexDat(const Index& i1, const QN& q1,
               const Args&... args) 
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
    fill() { }

    template<long J, typename... Args>
    void
    fill(const Index& i, const QN& q, const Args&... rest)
        {
#ifdef DEBUG
        assert(J < iq_.size());
#endif
        iq_[J] = IndexQN(i,q);
        fill<J+1>(rest...);
        }
    };

long
totalM(const IQIndexDat::storage& storage);

inline IQIndex::
IQIndex() 
    : 
    dir_(Neither)
    { }

namespace detail {

long inline
totalM()
    {
    return 0;
    }

template<typename... Args>
long
totalM(const Index& i, const QN& q,
       const Args&... args)
    {
    return i.m()+detail::totalM(args...);
    }

} //namespace detail


template<typename... Args>
IQIndex::
IQIndex(const std::string& name, 
        const Index& i1, const QN& q1, 
        const Args&... rest)
    : 
    Index(name,detail::totalM(i1,q1,rest...),i1.type(),i1.primeLevel()), 
    pd(std::make_shared<IQIndexDat>(i1,q1,rest...)),
    dir_(Out)
    { }

template<typename... Args>
IQIndex::
IQIndex(const std::string& name, 
        Arrow dir,
        const Index& i1, const QN& q1, 
        const Args&... rest)
    : 
    Index(name,detail::totalM(i1,q1,rest...),i1.type(),i1.primeLevel()), 
    pd(std::make_shared<IQIndexDat>(i1,q1,rest...)),
    dir_(dir)
    { }

inline IQIndex::
IQIndex(const std::string& name, 
        storage&& ind_qn, 
        Arrow dir, int plev) 
    : 
    Index(name,totalM(ind_qn),ind_qn.front().index.type(),plev),
    pd(std::make_shared<IQIndexDat>(std::move(ind_qn))),
    dir_(dir)
    { }


} //namespace itensor

#endif
