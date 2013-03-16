//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQINDEX_H
#define __ITENSOR_IQINDEX_H
#include "itensor.h"
#include "qn.h"

// Forward declarations
class IndexQN;
class IQIndexDat;
class IQIndexVal;

typedef boost::shared_ptr<IQIndexDat>
IQIndexDatPtr;

//
// IQIndex
//

class IQIndex : public Index
    {
    public:

    typedef std::vector<IndexQN>
    Storage;

    const Storage&
    indices() const;

    int 
    nindex() const;

    const Index& 
    index(int i) const;

    const QN& 
    qn(int i) const;

    int 
    primeLevel() const { return Index::primeLevel(); }
    void 
    primeLevel(int val);

    //------------------------------------------
    //IQIndex: Constructors

    IQIndex();

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
            std::vector<IndexQN>& ind_qn, 
            Arrow dir = Out, int plev = 0);

    IQIndex(const IQIndex& other, 
            std::vector<IndexQN>& ind_qn);

    IQIndex(const Index& other, 
            const Index& i1, const QN& q1, 
            Arrow dir = Out);


    explicit 
    IQIndex(const std::string& name,
            IndexType it = Link, 
            Arrow dir = Out, 
            int plev = 0);


    explicit 
    IQIndex(std::istream& s) { read(s); }

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

    //------------------------------------------
    //IQIndex: quantum number methods

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

    private:

    /////////////
    Arrow _dir;

    boost::shared_ptr<IQIndexDat> pd;
    /////////////

    IQIndex(const Index& index, const IQIndexDatPtr& pdat);

    void 
    solo();

    }; //class IQIndex


//
// IQIndexVal
//

class IQIndexVal : public IQIndex
    {
    public:

    int i;

    IQIndexVal();

    IQIndexVal(const IQIndex& iqindex, int i_);

    Index index() const;

    using IQIndex::index;

    QN qn() const;

    using IQIndex::qn;

    bool
    operator==(const IQIndexVal& other) const;

    operator IndexVal() const;

    IndexVal blockIndexVal() const;

    ITensor 
    operator*(const IndexVal& iv) const 
        { 
        return IndexVal(Index(*this),i) * iv; 
        }

    static const IQIndexVal& Null()
        {
        static const IQIndexVal Null_;
        return Null_;
        }

    private:

    void 
    calc_ind_ii(int& j, int& ii) const;

    };

//
// IndexQN
//

class IndexQN : public Index
    {
    public:

    QN qn;

    IndexQN() { }

    IndexQN(const Index& i, const QN& q) : Index(i), qn(q) { }

    void 
    write(std::ostream& s) const { Index::write(s); qn.write(s); }
    void 
    read(std::istream& s) { Index::read(s); qn.read(s); }
    };


std::string 
showm(const IQIndex& I);

int 
offset(const IQIndex& I, const Index& i);


std::ostream& 
operator<<(std::ostream &o, const IQIndex &I);

std::ostream& 
operator<<(std::ostream &s, const IndexQN& x);

std::ostream& 
operator<<(std::ostream& s, const IQIndexVal& iv);

#endif
