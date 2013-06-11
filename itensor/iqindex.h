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
    //
    //Constructors
    //

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
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            const Index& i4, const QN& q4,
            const Index& i5, const QN& q5,
            const Index& i6, const QN& q6,
            Arrow dir = Out);

    IQIndex(const std::string& name, 
            std::vector<IndexQN>& ind_qn, 
            Arrow dir = Out, 
            int plev = 0);

    //
    //Accessor Methods
    //

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

    Arrow 
    dir() const { return dir_; }

    int 
    primeLevel() const { return Index::primeLevel(); }
    void 
    primeLevel(int val);

    //
    // Prime level methods
    //

    void 
    prime(int inc = 1);

    void 
    prime(IndexType type, int inc = 1);

    void 
    noprime(IndexType type = All);

    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //
    // Operators
    //

    IQIndexVal 
    operator()(int n) const;

    //
    // Other methods
    //

    void 
    conj();

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    //
    // Static Index instances
    //

    static const 
    IQIndex& Null();

    private:

    /////////////
    Arrow dir_;

    boost::shared_ptr<IQIndexDat> pd;
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

class IQIndexVal : public IQIndex
    {
    public:

    int i;

    IQIndexVal();

    IQIndexVal(const IQIndex& iqindex, int i_);

    IndexQN
    indexqn() const;

    const QN&
    qn() const;

    const QN&
    qn(int j) const { return IQIndex::qn(j); }

    bool
    operator==(const IQIndexVal& other) const;

    operator IndexVal() const;

    IndexVal 
    blockIndexVal() const;

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

    };


bool
hasindex(const IQIndex& I, const Index& i);

int 
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

#endif
