//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __IQINDEX_H
#define __IQINDEX_H
#include "index.h"
#include "itensor.h"

/*
* Conventions regarding arrows:
*
* * Arrows point In or Out, never right/left/up/down.
*
* * The Site indices of a ket point Out.
*
* * Conjugation switches arrow directions.
*
* * All arrows flow Out from the ortho center of an MPS 
*   (assuming it's a ket - In if it's a bra).
*
* * IQMPOs are created with the same arrow structure as if they are 
*   orthogonalized to site 1, but this is just a default since they 
*   aren't actually ortho. If position is called on an IQMPO it follows 
*   the same convention as for an MPS except Site indices point In and 
*   Site' indices point Out.
*
* * Local site operators have two IQIndices, one unprimed and pointing In, 
*   the other primed and pointing Out.
*
*/

// Forward declarations
class QN;
struct inqn;
class IQIndexDat;
struct IQIndexVal;


//
// IQIndex
//

class IQIndex : public Index
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

    IQIndex(PrimeType pt, const IQIndex& other, int inc = 1);

    explicit 
    IQIndex(std::istream& s);

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    IQIndex(Imaker im);

    static const IQIndex& 
    Null()
        {
        static const IQIndex Null_(makeNull);
        return Null_;
        }

    static const IQIndex& 
    IndReIm()
        {
        static const IQIndex IndReIm_(makeReIm);
        return IndReIm_;
        }

    static const IQIndex& 
    IndReImP()
        {
        static const IQIndex IndReImP_(makeReImP);
        return IndReImP_;
        }

    static const IQIndex& 
    IndReImPP()
        {
        static const IQIndex IndReImPP_(makeReImPP);
        return IndReImPP_;
        }

    IQIndexVal 
    operator()(int n) const;

    //------------------------------------------
    //IQIndex: methods for querying m's

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
    doprime(PrimeType pt, int inc = 1);

    void 
    mapprime(int plevold, int plevnew, PrimeType pt = primeBoth);

    void 
    noprime(PrimeType pt = primeBoth);

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
        { return IQIndex(primeBoth,I,inc); }

    friend std::ostream& 
    operator<<(std::ostream &o, const IQIndex &I);

    void 
    print(std::string name = "") const;

    typedef std::vector<inqn>::iterator 
    iq_it;
    typedef std::vector<inqn>::const_iterator 
    const_iq_it;

    private:

    Arrow _dir;

    boost::intrusive_ptr<IQIndexDat> pd;

    void 
    solo();


    }; //class IQIndex




//
// QN
//

class QN
    {
    public:

    QN(int sz=0,int Nf=0) 
        : _sz(sz), _Nf(Nf), _Nfp(abs(Nf%2)) 
        { }

    QN(int sz,int Nf,int Nfp) 
        : _sz(sz), _Nf(Nf), _Nfp(abs(Nfp%2))
        { assert(_Nf==0 || abs(_Nf%2) == _Nfp); }

    QN(std::istream& s) { read(s); }

    int 
    sz() const { return _sz; }

    int 
    Nf() const { return _Nf; }

    int 
    Nfp() const { assert(_Nfp == 0 || _Nfp == 1); return _Nfp; }

    int 
    sign() const { return (_Nfp == 0 ? +1 : -1); }

    void 
    write(std::ostream& s) const 
        { 
        s.write((char*)&_sz,sizeof(_sz)); 
        s.write((char*)&_Nf,sizeof(_Nf)); 
        s.write((char*)&_Nfp,sizeof(_Nfp)); 
        }

    void 
    read(std::istream& s) 
        { 
        s.read((char*)&_sz,sizeof(_sz)); 
        s.read((char*)&_Nf,sizeof(_Nf)); 
        s.read((char*)&_Nfp,sizeof(_Nfp)); 
        }

    QN 
    operator+(const QN &other) const
        { QN res(*this); res+=other; return res; }

    QN 
    operator-(const QN &other) const
        { QN res(*this); res-=other; return res; }

    QN& 
    operator+=(const QN &other)
        {
        _sz+=other._sz; _Nf+=other._Nf; _Nfp = abs(_Nfp+other._Nfp)%2;
        return *this;
        }

    QN& 
    operator-=(const QN &other)
        {
        _sz-=other._sz; _Nf-=other._Nf; _Nfp = abs(_Nfp-other._Nfp)%2;
        return *this;
        }

    QN 
    operator-() const  
        { return QN(-_sz,-_Nf,_Nfp); }
    
    QN 
    negated() const { return QN(-_sz,-_Nf,_Nfp); }

    //Multiplication and division should only be used to change the sign
    QN& 
    operator*=(int i) { assert(i*i == 1); _sz*=i; _Nf*=i; return *this; }

    QN 
    operator*(int i) const { QN res(*this); res*=i; return res; }

    QN 
    operator/(int i) const { QN res(*this); res*=i; return res; }

    std::string 
    toString() const
        { return (boost::format("(%+d:%d:%s)")%_sz%_Nf%(_Nfp==1 ? "-" : "+")).str(); }

    inline friend std::ostream& 
    operator<<(std::ostream &o, const QN &q)
        { return o<< boost::format("sz = %d, Nf = %d, fp = %s") % q.sz() % q.Nf() % (q.sign() < 0 ? "-" : "+"); }

    void 
    print(std::string name = "") const
        { std::cerr << "\n" << name << " =\n" << *this << "\n"; }

    private:

    int _sz, 
        _Nf, 
        _Nfp; //_Nfp stands for fermion number parity, and tracks whether Nf is even or odd
    };

inline bool 
operator==(const QN &a,const QN &b)
    { return a.sz() == b.sz() && a.Nf() == b.Nf() && a.Nfp() == b.Nfp(); }

inline bool 
operator!=(const QN &a,const QN &b)
    { return a.sz() != b.sz() || a.Nf() != b.Nf() || a.Nfp() != b.Nfp(); }

inline bool 
operator<(const QN &a,const QN &b)
    { return a.sz() < b.sz() || (a.sz() == b.sz() && a.Nf() < b.Nf()) 
             || (a.sz() == b.sz() && a.Nf() == b.Nf() && a.Nfp() < b.Nfp()); }

inline QN 
operator*(int i,const QN& a)
    { return a*i; }




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

    std::vector<inqn> iq_;

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
    IQIndexDat(const IQIndexDat& other);

    explicit 
    IQIndexDat(std::istream& s);

    explicit 
    IQIndexDat(Imaker im);

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    static IQIndexDat* Null()
        {
        static IQIndexDat Null_(makeNull);
        return &Null_;
        }

    static IQIndexDat* ReImDat()
        {
        static IQIndexDat ReImDat_(makeReIm);
        return &ReImDat_;
        }

    static IQIndexDat* ReImDatP()
        {
        static IQIndexDat ReImDatP_(makeReImP);
        return &ReImDatP_;
        }

    static IQIndexDat* ReImDatPP()
        {
        static IQIndexDat ReImDatPP_(makeReImPP);
        return &ReImDatPP_;
        }

    friend void intrusive_ptr_add_ref(IQIndexDat* p);
    friend void intrusive_ptr_release(IQIndexDat* p);

    int 
    count() const { return numref; }

    typedef std::vector<inqn>::iterator 
    iq_it;
    typedef std::vector<inqn>::const_iterator 
    const_iq_it;

    private:

    //Disallow copying
    void 
    operator=(const IQIndexDat&);

    mutable unsigned int numref;

    const bool is_static_;
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

    IQIndexVal primed() const 
        { return IQIndexVal(iqind.primed(),i); }

    void conj() 
        { iqind.conj(); }

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
        static const IQIndexVal Null_(IQIndex::Null(),1);
        return Null_;
        }

    private:

    void 
    calc_ind_ii(int& j, int& ii) const;

    };


class DoPrimer // Functor which applies doprime within STL's for_each, etc
    {
    public:

    PrimeType pt; 

    int inc;

    DoPrimer (PrimeType _pt, int _inc = 1) 
        : pt(_pt), 
          inc(_inc) 
        { }

    void 
    operator()(inqn& iq) const { iq.index.doprime(pt,inc); }
    void 
    operator()(Index& i) const { i.doprime(pt,inc); }
    void 
    operator()(ITensor& it) const { it.doprime(pt,inc); }
    void 
    operator()(IQIndex &iqi) const { iqi.doprime(pt,inc); }
    };

class MapPrimer // Functor which applies mapprime within STL's for_each, etc
    {
    public:

    PrimeType pt;

    int plevold, plevnew;

    MapPrimer (int _plevold,int _plevnew,PrimeType _pt = primeBoth) 
		: pt(_pt), 
          plevold(_plevold), 
          plevnew(_plevnew) 
        {}

    void 
    operator()(inqn& iq) const { iq.index.mapprime(plevold,plevnew,pt); }
    void 
    operator()(Index& i) const { i.mapprime(plevold,plevnew,pt); }
    void 
    operator()(ITensor& it) const { it.mapprime(plevold,plevnew,pt); }
    void 
    operator()(IQIndex &iqi) const { iqi.mapprime(plevold,plevnew,pt); }
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
