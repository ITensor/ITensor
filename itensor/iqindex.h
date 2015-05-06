//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQINDEX_H
#define __ITENSOR_IQINDEX_H
#include "index.h"
#include "qn.h"

namespace itensor {

// Forward declarations
class ITensor;
class IndexQN;
class IQIndexDat;
class IQIndexVal;

using IQIndexDatPtr = shared_ptr<IQIndexDat>;

//
// IQIndex
//

class IQIndex : public Index
    {
    public:

    using storage = std::vector<IndexQN>;
    //
    //Constructors
    //

    IQIndex();

    template<typename... Rest>
    IQIndex(const std::string& name, 
            const Index& i1, const QN& q1, 
            const Rest&... etc);

    template<typename... Rest>
    IQIndex(const std::string& name, 
            Arrow dir,
            const Index& i1, const QN& q1, 
            const Rest&... etc);

    IQIndex(const std::string& name, 
            storage&& ind_qn, 
            Arrow dir = Out, 
            int plev = 0);

    //
    //Accessor Methods
    //


    const storage&
    inds() const;

    long 
    nindex() const;

    //1-indexed
    const Index& 
    index(long i) const;

    //0-indexed
    const Index& 
    operator[](long i) const;

    //1-indexed
    const QN& 
    qn(long i) const;

    Arrow 
    dir() const { return dir_; }
    //void 
    //dir(Arrow ndir) { dir_ = ndir; }

    int 
    primeLevel() const { return Index::primeLevel(); }
    IQIndex& 
    primeLevel(int val);

    //
    // Prime level methods
    //

    IQIndex& 
    prime(int inc = 1);

    IQIndex& 
    prime(IndexType type, int inc = 1);

    IQIndex& 
    noprime(IndexType type = All);

    IQIndex& 
    mapprime(int plevold, int plevnew, IndexType type = All);

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

    void 
    write(std::ostream& s) const;

    IQIndex& 
    read(std::istream& s);

    private:

    /////////////
    IQIndexDatPtr pd;
    Arrow dir_;
    /////////////

    IQIndex(const Index& index, const IQIndexDatPtr& pdat);

    void 
    solo();

    }; //class IQIndex

//
// IndexQN
//

class IndexQN : public Index
    {
    public:
    
    using IndexValT = IQIndexVal;

    QN qn;

    IndexQN() { }

    IndexQN(const Index& i, const QN& q) : Index(i), qn(q) { }

    void 
    write(std::ostream& s) const { Index::write(s); qn.write(s); }
    void 
    read(std::istream& s) { Index::read(s); qn.read(s); }
    };


//
// IQIndexVal
//

class IQIndexVal
    {
    public:

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

    static const IQIndexVal& Null()
        {
        static const IQIndexVal Null_;
        return Null_;
        }

    };


IQIndex inline
dag(IQIndex res) { res.dag(); return res; }

IndexQN inline
dag(IndexQN res) { res.dag(); return res; }

IQIndexVal inline
dag(IQIndexVal res) { res.dag(); return res; }

IQIndex inline
operator^(IQIndex I, int inc) { I.prime(inc); return I; }

IndexQN inline
operator^(IndexQN I, int inc) { I.prime(inc); return I; }

IQIndexVal inline
operator^(IQIndexVal I, int inc) { I.prime(inc); return I; }

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
    IQIndexDat(storage&& ind_qn) { iq_ = std::move(ind_qn); }

    const storage&
    inds() const { return iq_; }

    long
    size() { return iq_.size(); }

    const Index&
    index(long i) { return iq_[i-1]; }

    const Index&
    operator[](long i) { return iq_[i]; }

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
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    static const IQIndexDatPtr& Null();

    void
    makeCopyOf(const IQIndexDat& other) { iq_ = other.iq_; }

    private:

    //////////////////

    storage iq_;

    /////////////////

    //Disallow copying using =
    void 
    operator=(const IQIndexDat&);

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

}; //namespace detail


template<typename... Args>
IQIndex::
IQIndex(const std::string& name, 
        const Index& i1, const QN& q1, 
        const Args&... rest)
    : 
    Index(name,detail::totalM(i1,q1,rest...),i1.type(),i1.primeLevel()), 
    pd(make_shared<IQIndexDat>(i1,q1,rest...)),
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
    pd(make_shared<IQIndexDat>(i1,q1,rest...)),
    dir_(dir)
    { }

inline IQIndex::
IQIndex(const std::string& name, 
        storage&& ind_qn, 
        Arrow dir, int plev) 
    : 
    Index(name,totalM(ind_qn),ind_qn.front().type(),plev),
    pd(make_shared<IQIndexDat>(std::move(ind_qn))),
    dir_(dir)
    { 
    }


}; //namespace itensor

#endif
