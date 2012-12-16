//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __IQINDEX_H
#define __IQINDEX_H
#include "index.h"
#include "itensor.h"
#include "qn.h"


// Forward declarations
struct inqn;
class IQIndexDat;
struct IQIndexVal;



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
    operator<(const IQIndex& other) const
        { return index_.operator<(other.index_); }

    bool 
    noprime_equals(const IQIndex& other) const
        { return index_.noprime_equals(other.index_); }

    //------------------------------------------
    //IQIndex: methods for querying m's

    int
    m() const { return index_.m(); }

    int 
    biggestm() const;

    std::string 
    showm() const;

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
    prime(int inc = 1) { prime(All,inc); }

    void 
    prime(IndexType type, int inc = 1);

    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    void 
    noprime(IndexType type = All);

    IQIndex friend inline
    noprime(const IQIndex& I)
        { 
        IQIndex J(I);
        J.noprime(); 
        return J; 
        }

    IQIndex 
    primed(int inc = 1) const;

    friend inline IQIndex
    primed(const IQIndex& I, int inc = 1)
        { return IQIndex(All,I,inc); }

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

    boost::intrusive_ptr<IQIndexDat> pd;

    explicit
    IQIndex(Index::Imaker im);

    void 
    solo();


    }; //class IQIndex





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



//
// IQIndexDat
//

class IQIndexDat
    {
    public:


    IQIndexDat();

    IQIndexDat(const Index& i1, const QN& q1);

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2);

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3);

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3,
               const Index& i4, const QN& q4);

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3,
               const Index& i4, const QN& q4,
               const Index& i5, const QN& q5);

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3,
               const Index& i4, const QN& q4,
               const Index& i5, const QN& q5,
               const Index& i6, const QN& q6);

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3,
               const Index& i4, const QN& q4,
               const Index& i5, const QN& q5,
               const Index& i6, const QN& q6,
               const Index& i7, const QN& q7);

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3,
               const Index& i4, const QN& q4,
               const Index& i5, const QN& q5,
               const Index& i6, const QN& q6,
               const Index& i7, const QN& q7,
               const Index& i8, const QN& q8);

    IQIndexDat(std::vector<inqn>& ind_qn);

    explicit 
    IQIndexDat(std::istream& s);

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    static IQIndexDat* Null();

    static IQIndexDat* ReImDat();

    static IQIndexDat* ReImDatP();

    static IQIndexDat* ReImDatPP();

    friend void 
    intrusive_ptr_add_ref(IQIndexDat* p);

    friend void 
    intrusive_ptr_release(IQIndexDat* p);

    int 
    count() const { return numref; }

    typedef std::vector<inqn>::iterator 
    iq_it;

    typedef std::vector<inqn>::const_iterator 
    const_iq_it;

    private:

    friend class IQIndex;

    explicit 
    IQIndexDat(Index::Imaker im);

    //////////////////
    //
    // Data Members
    //

    std::vector<inqn> iq_;

    mutable int numref;

    const bool is_static_;

    //
    /////////////////

    explicit 
    IQIndexDat(const IQIndexDat& other);

    //Disallow copying using =
    void 
    operator=(const IQIndexDat&);

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

    friend inline IQIndexVal
    primed(const IQIndexVal& IV, int inc = 1)
        {
        return IQIndexVal(primed(IV.iqind,inc),IV.i); 
        }

    void conj() 
        { iqind.conj(); }

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


class Primer // Functor which applies prime within STL's for_each, etc
    {
    public:

    IndexType type; 

    int inc;

    Primer (IndexType type_, int inc_ = 1) 
        : type(type_), 
          inc(inc_) 
        { }

    void 
    operator()(inqn& iq) const { iq.index.prime(type,inc); }
    void 
    operator()(Index& i) const { i.prime(type,inc); }
    void 
    operator()(ITensor& it) const { it.prime(type,inc); }
    void 
    operator()(IQIndex &iqi) const { iqi.prime(type,inc); }
    };

class MapPrimer // Functor which applies mapprime within STL's for_each, etc
    {
    public:

    IndexType type;

    int plevold, plevnew;

    MapPrimer (int _plevold,int _plevnew,IndexType _type = All) 
		: type(_type), 
          plevold(_plevold), 
          plevnew(_plevnew) 
        {}

    void 
    operator()(inqn& iq) const { iq.index.mapprime(plevold,plevnew,type); }
    void 
    operator()(Index& i) const { i.mapprime(plevold,plevnew,type); }
    void 
    operator()(ITensor& it) const { it.mapprime(plevold,plevnew,type); }
    void 
    operator()(IQIndex &iqi) const { iqi.mapprime(plevold,plevnew,type); }
    };

class IndEq // Functor which checks if the index is equal to a specified value within STL's for_each, etc
    {
    public:

    Index i;

    IndEq(Index _i) 
        : i(_i) {}

    bool 
    operator()(const inqn &j) const { return i == j.index; }
    };


#endif
