//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __IQINDEX_H
#define __IQINDEX_H
#include "itensor.h"
#include "qn.h"


// Forward declarations
struct inqn;
class IQIndexDat;
struct IQIndexVal;

typedef boost::shared_ptr<IQIndexDat>
IQIndexDatPtr;


//
// IQIndex
//

class IQIndex
    {
    public:

    const std::vector<inqn>& 
    iq() const;

    int 
    nindex() const;

    const Index& 
    index(int i) const;

    const QN& 
    qn(int i) const;

    IndexType 
    type() const { return index_.type(); }

    std::string 
    name() const { return index_.name(); }

    const std::string&
    rawname() const { return index_.rawname(); }

    Real 
    uniqueReal() const { return index_.uniqueReal(); }

    int 
    primeLevel() const { return index_.primeLevel(); }
    void 
    primeLevel(int val) { index_.primeLevel(val); }

    bool
    isNull() const { return index_.isNull(); }
    bool
    isNotNull() const { return index_.isNotNull(); }

    //------------------------------------------
    //IQIndex: Constructors

    IQIndex();

    explicit 
    IQIndex(const Index& other, Arrow dir = Out);

    explicit 
    IQIndex(const std::string& name,
            IndexType it = Link, 
            Arrow dir = Out, 
            int plev = 0);

    IQIndex(const std::string& name, 
            const Index& i1, const QN& q1, 
            Arrow dir = Out);

    IQIndex(const std::string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            Arrow dir = Out);

    IQIndex(const std::string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            Arrow dir = Out);

    IQIndex(const std::string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            const Index& i4, const QN& q4,
            Arrow dir = Out);

    IQIndex(const std::string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            const Index& i4, const QN& q4,
            const Index& i5, const QN& q5,
            Arrow dir = Out);

    IQIndex(const std::string& name, 
            std::vector<inqn>& ind_qn, 
            Arrow dir = Out, int plev = 0);

    IQIndex(const IQIndex& other, 
            std::vector<inqn>& ind_qn);

    IQIndex(const Index& other, 
            const Index& i1, const QN& q1, 
            Arrow dir = Out);

    IQIndex(IndexType type, const IQIndex& other, int inc = 1);

    explicit 
    IQIndex(std::istream& s);

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    static const IQIndex& Null();

    static const IQIndex& IndReIm();

    static const IQIndex& IndReImP();

    static const IQIndex& IndReImPP();

    //------------------------------------------
    //IQIndex: operators

    IQIndexVal 
    operator()(int n) const;

    operator Index() const { return index_; }

    bool 
    operator==(const IQIndex& other) const
        { return index_.operator==(other.index_); }

    bool 
    operator!=(const IQIndex& other) const
        { return index_.operator!=(other.index_); }

    bool 
    operator<(const IQIndex& other) const
        { return index_.operator<(other.index_); }

    bool 
    noprimeEquals(const IQIndex& other) const
        { return index_.noprimeEquals(other.index_); }

    //------------------------------------------
    //IQIndex: methods for querying m's

    int
    m() const { return index_.m(); }

    int 
    biggestm() const;

    //------------------------------------------
    //IQIndex: quantum number methods

    void 
    negate();

    friend IQIndex 
    negate(IQIndex I); // Quantum numbers negated
     
    QN 
    qn(const Index& i) const;

    Arrow 
    dir() const;

    void 
    conj();

    //------------------------------------------
    //IQIndex: index container methods


    const Index& 
    findbyqn(QN q) const;

    bool 
    hasindex(const Index& i) const;

    bool 
    hasindex_noprime(const Index& i) const;

    int 
    offset(const Index& I) const;

    //------------------------------------------
    //IQIndex: prime methods

    void 
    prime(int inc = 1);

    void 
    prime(IndexType type, int inc = 1);

    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    void 
    noprime(IndexType type = All);

    void 
    prime(int inc = 1) const;

    friend std::ostream& 
    operator<<(std::ostream &o, const IQIndex &I);

    void 
    print(std::string name = "") const;

    typedef std::vector<inqn>::iterator 
    iq_it;
    typedef std::vector<inqn>::const_iterator 
    const_iq_it;

    private:

    Index index_;

    Arrow _dir;

    boost::shared_ptr<IQIndexDat> pd;

    explicit
    IQIndex(Index::Imaker im);

    void 
    solo();


    }; //class IQIndex


IQIndex inline
primed(IQIndex I, int inc = 1) { I.prime(inc); return I; }

IQIndex inline
primed(IQIndex I, IndexType type, int inc = 1) { I.prime(type,inc); return I; }

// Return a copy of this Index with primelevel set to zero.
IQIndex inline
deprimed(IQIndex I) { I.noprime();  return I; }

std::string 
showm(const IQIndex& I);


//
// inqn
//

struct inqn
    {
    Index index;
    QN qn;
    inqn() { }
    inqn(const Index& i, QN q) 
        : index(i), 
          qn(q) 
        { }

    void 
    write(std::ostream& s) const { index.write(s); qn.write(s); }

    void 
    read(std::istream& s) { index.read(s); qn.read(s); }

    inline friend std::ostream& 
    operator<<(std::ostream &o, const inqn& x)
        { o << "inqn: " << x.index << " (" << x.qn << ")\n"; return o; }
    };


enum IQmaker {makeSing};



//
// IQIndexVal
//

struct IQIndexVal
    {
    IQIndex iqind; 

    int i;

    IQIndexVal();

    IQIndexVal(const IQIndex& iqindex, int i_);

    Index index() const;

    QN qn() const;

    void
    prime(int inc = 1) { iqind.prime(inc); }

    void
    prime(IndexType type, int inc = 1) { iqind.prime(type,inc); }

    void 
    conj() { iqind.conj(); }

    bool
    operator==(const IQIndexVal& other) const;

    operator IndexVal() const;

    IndexVal blockIndexVal() const;

    //operator ITensor() const;

    ITensor 
    operator*(const IndexVal& iv) const 
        { 
        return IndexVal(Index(iqind),i) * iv; 
        }

    void 
    print(std::string name = "") const
        { std::cerr << "\n" << name << " =\n" << *this << "\n"; }

    inline friend std::ostream& 
    operator<<(std::ostream& s, const IQIndexVal& iv)
        { return s << "IQIndexVal: i = " << iv.i << ", iqind = " << iv.iqind << "\n"; }

    static const IQIndexVal& Null()
        {
        static const IQIndexVal Null_(Index::makeNull);
        return Null_;
        }

    private:

    void 
    calc_ind_ii(int& j, int& ii) const;

    explicit
    IQIndexVal(Index::Imaker im);

    };




#endif
