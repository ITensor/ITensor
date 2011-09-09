#ifndef __IQINDEX_H
#define __IQINDEX_H
#include "index.h"

/*
* Conventions regarding arrows:
*
* * Arrows point In or Out, never right/left/up/down.
*
* * The Site indices of a ket point out.
*
* * Conjugation switches arrow directions.
*
* * All arrows flow out from the ortho center of an MPS (ket - in if it's a bra).
*
* * IQMPOs are created with the same arrow structure as if they are orthogonalized
*   to site 1, but this is just a default since they aren't actually ortho. If position 
*   is called on an IQMPO it follows the same convention as for an MPS except Site 
*   indices point In and Site' indices point Out.
*
* * The Virtual IQIndex moves with the orthogonality center. If two IQTensors are
*   multiplied, the QNs of their Virtual indices are added together if they both have one.
*
*/

class QN
{
    int _sz, _Nf, _Nfp; //_Nfp stands for fermion number parity, and tracks whether Nf is even or odd
public:
    QN(int sz=0,int Nf=0) : _sz(sz), _Nf(Nf), _Nfp(abs(Nf%2)) { }
    QN(int sz,int Nf,int Nfp) : _sz(sz), _Nf(Nf), _Nfp(abs(Nfp%2))
    { assert(_Nf==0 || abs(_Nf%2) == _Nfp); }

    QN(istream& s) { read(s); }

    int sz() const { return _sz; }
    int Nf() const { return _Nf; }
    int Nfp() const { assert(_Nfp == 0 || _Nfp == 1); return _Nfp; }
    int fp() const { return (_Nfp == 0 ? +1 : -1); }
    //int& sz() { return _sz; }
    //int& Nf() { return _Nf; }
    //int& Nfp() { assert(_Nfp == 0 || _Nfp == 1); return _Nfp; }

    void write(ostream& s) const 
    { 
        s.write((char*)&_sz,sizeof(_sz)); 
        s.write((char*)&_Nf,sizeof(_Nf)); 
        s.write((char*)&_Nfp,sizeof(_Nfp)); 
    }
    void read(istream& s) 
    { 
        s.read((char*)&_sz,sizeof(_sz)); 
        s.read((char*)&_Nf,sizeof(_Nf)); 
        s.read((char*)&_Nfp,sizeof(_Nfp)); 
    }

    QN operator+(const QN &other) const
	{ QN res(*this); res+=other; return res; }
    QN operator-(const QN &other) const
	{ QN res(*this); res-=other; return res; }
    QN& operator+=(const QN &other)
	{
        _sz+=other._sz; _Nf+=other._Nf; _Nfp = abs(_Nfp+other._Nfp)%2;
        return *this;
	}
    QN& operator-=(const QN &other)
	{
        _sz-=other._sz; _Nf-=other._Nf; _Nfp = abs(_Nfp-other._Nfp)%2;
        return *this;
	}
    QN operator-() const  
	{ return QN(-_sz,-_Nf,_Nfp); }
    
    QN negated() const { return QN(-_sz,-_Nf,_Nfp); }

    //Multiplication and division should only be used to change the sign
    QN& operator*=(int i) { assert(i*i == 1); _sz*=i; _Nf*=i; return *this; }
    QN operator*(int i) const { QN res(*this); res*=i; return res; }
    QN operator/(int i) const { QN res(*this); res*=i; return res; }

    inline string toString() const
    {  return (format("(%+d:%d:%s)")%_sz%_Nf%(_Nfp==1 ? "-" : "+")).str(); }

    inline friend ostream& operator<<(ostream &o, const QN &q)
    { return o<< format("sz = %d, Nf = %d, fp = %s") % q.sz() % q.Nf() % (q.fp() < 0 ? "-" : "+"); }

    void print(string name = "") const
    { cerr << "\n" << name << " =\n" << *this << "\n"; }
};
inline bool operator==(const QN &a,const QN &b)
    { return a.sz() == b.sz() && a.Nf() == b.Nf() && a.Nfp() == b.Nfp(); }
inline bool operator!=(const QN &a,const QN &b)
    { return a.sz() != b.sz() || a.Nf() != b.Nf() || a.Nfp() != b.Nfp(); }
inline bool operator<(const QN &a,const QN &b)
    { return a.sz() < b.sz() || (a.sz() == b.sz() && a.Nf() < b.Nf()) 
             || (a.sz() == b.sz() && a.Nf() == b.Nf() && a.Nfp() < b.Nfp()); }
inline QN operator*(int i,const QN& a)
    { return a*i; }


struct inqn
{
    Index index;
    QN qn;
    inqn() { }
    inqn(const Index& i, QN q) : index(i), qn(q) { }

    void write(ostream& s) const { index.write(s); qn.write(s); }
    void read(istream& s) { index.read(s); qn.read(s); }
    inline friend ostream& operator<<(ostream &o, const inqn& x)
    { o << "inqn: " << x.index << " (" << x.qn << ")\n"; return o; }
};

class IQIndex;
class DoPrimer		// Functor which applies doprime within STL's for_each, etc
{
public:
    PrimeType pt; int inc;
    DoPrimer (PrimeType _pt, int _inc = 1) : pt(_pt), inc(_inc) {}
    void operator()(inqn& iq) const { iq.index.doprime(pt,inc); }
    void operator()(Index& i) const { i.doprime(pt,inc); }
    void operator()(ITensor& it) const { it.doprime(pt,inc); }
    void operator()(IQIndex& iqi) const;
};
class MapPrimer // Functor which applies mapprime within STL's for_each, etc
{
public:
    PrimeType pt;
    int plevold, plevnew;
    MapPrimer (int _plevold,int _plevnew,PrimeType _pt = primeBoth) 
		: pt(_pt), plevold(_plevold), plevnew(_plevnew) {}
    void operator()(inqn& iq) const { iq.index.mapprime(plevold,plevnew,pt); }
    void operator()(Index& i) const { i.mapprime(plevold,plevnew,pt); }
    void operator()(ITensor& it) const { it.mapprime(plevold,plevnew,pt); }
    void operator()(IQIndex& iqi) const;
};
class IndEq // Functor which checks if the index is equal to a specified value within STL's for_each, etc
{
public:
    Index i;
    IndEq(Index _i) : i(_i) {}
    bool operator()(const inqn &j) const { return i == j.index; }
};

class IQIndexDat
{
    mutable unsigned int numref;
    const bool is_static_;
public:
    vector<inqn> iq_;

    IQIndexDat() : numref(0), is_static_(false) { }

    IQIndexDat(const Index& i1, const QN& q1)
    : numref(0), is_static_(false)
    {
        iq_.push_back(inqn(i1,q1));
    }

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2)
    : numref(0), is_static_(false)
    {
        iq_.push_back(inqn(i1,q1));
        iq_.push_back(inqn(i2,q2));
    }

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3)
    : numref(0), is_static_(false)
    {
        iq_.push_back(inqn(i1,q1));
        iq_.push_back(inqn(i2,q2));
        iq_.push_back(inqn(i3,q3));
    }

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2, const QN& q2,
               const Index& i3, const QN& q3,
               const Index& i4, const QN& q4)
    : numref(0), is_static_(false)
    {
        iq_.push_back(inqn(i1,q1));
        iq_.push_back(inqn(i2,q2));
        iq_.push_back(inqn(i3,q3));
        iq_.push_back(inqn(i4,q4));
    }

    IQIndexDat(vector<inqn>& ind_qn)
    : numref(0), is_static_(false)
    { iq_.swap(ind_qn); }

    explicit IQIndexDat(const IQIndexDat& other) : numref(0), is_static_(false), iq_(other.iq_)
    { }

    explicit IQIndexDat(istream& s) : numref(0), is_static_(false) { read(s); }

    void write(ostream& s) const
    {
        size_t size = iq_.size();
        s.write((char*)&size,sizeof(size));
        foreach(const inqn& x,iq_) x.write(s);
    }

    void read(istream& s)
    {
        size_t size; s.read((char*)&size,sizeof(size));
        iq_.resize(size);
        foreach(inqn& x,iq_) x.read(s);
    }

    explicit IQIndexDat(Imaker im)
    : numref(0), is_static_(true)
    { 
        iq_.push_back(inqn(Index(im),QN())); 
    }

    friend inline void intrusive_ptr_add_ref(IQIndexDat* p) { ++(p->numref); }
    friend inline void intrusive_ptr_release(IQIndexDat* p) { if(!p->is_static_ && --(p->numref) == 0){ delete p; } }
    int count() const { return numref; }
private:
    void operator=(const IQIndexDat&);
};

extern IQIndexDat IQIndexDatNull, IQIndReDat, IQIndReDatP, IQIndReDatPP, IQIndEmptyVDat;

struct IQIndexVal;

class IQIndex : public Index
{
    Arrow _dir;
    intrusive_ptr<IQIndexDat> pd;
    void solo()
    {
        assert(pd != 0);
        if(pd->count() != 1)
        {
            intrusive_ptr<IQIndexDat> new_pd(new IQIndexDat(*pd));
            //new_pd->iq_ = pd->iq_;
            pd.swap(new_pd);
        }
    }
public:
    const vector<inqn>& iq() const { assert(pd != 0); return pd->iq_; }
    int nindex() const { return (int) pd->iq_.size(); }
    const Index& index(int i) const { assert(pd != 0); return GET(pd->iq_,i-1).index; }
    const QN& qn(int i) const { assert(pd != 0); return GET(pd->iq_,i-1).qn; }

    //------------------------------------------
    //IQIndex: Constructors

    IQIndex() : _dir(Out), pd(0) {}

    explicit IQIndex(const Index& other, Arrow dir = Out) : Index(other), _dir(dir), pd(0) {}

    explicit IQIndex(const string& name,IndexType it = Link, Arrow dir = Out, int plev = 0) : Index(name,1,it,plev), _dir(dir), pd(0) {}

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            Arrow dir = Out) 
    : Index(name,i1.m(),i1.type()), _dir(dir), pd(new IQIndexDat(i1,q1))
    {
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            Arrow dir = Out) 
    : Index(name,i1.m()+i2.m(),i1.type()), _dir(dir), 
    pd(new IQIndexDat(i1,q1,i2,q2))
    {
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            Arrow dir = Out) 
    : Index(name,i1.m()+i2.m()+i3.m(),i1.type()), _dir(dir),
    pd(new IQIndexDat(i1,q1,i2,q2,i3,q3))
    {
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            const Index& i4, const QN& q4,
            Arrow dir = Out) 
    : Index(name,i1.m()+i2.m()+i3.m()+i4.m(),i1.type()), _dir(dir),
    pd(new IQIndexDat(i1,q1,i2,q2,i3,q3,i4,q4))
    {
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, vector<inqn>& ind_qn, Arrow dir = Out, int plev = 0) 
    : Index(name,0,ind_qn.back().index.type(),plev), _dir(dir),
    pd(new IQIndexDat(ind_qn))
    { 
        //int* pm = const_cast<int*>(&(p->m_));
        //foreach(const inqn& x, pd->iq_) *pm += x.index.m();
        int mm = 0;
        foreach(const inqn& x, pd->iq_) mm += x.index.m();
        set_m(mm);
        setPrimeLevel(pd->iq_.back().index.primelevel);
    }

    IQIndex(const IQIndex& other, vector<inqn>& ind_qn)
    : Index(other.name(),0,other.type()), _dir(other._dir),
    pd(new IQIndexDat(ind_qn))
    { 
        //int* pm = const_cast<int*>(&(p->m_));
        //foreach(const inqn& x, pd->iq_) *pm += x.index.m();
        int mm = 0;
        foreach(const inqn& x, pd->iq_) mm += x.index.m();
        set_m(mm);
        setPrimeLevel(pd->iq_.back().index.primelevel);
    }

    IQIndex(const Index& other, 
            const Index& i1, const QN& q1, 
            Arrow dir = Out) 
    : Index(other), _dir(dir), pd(new IQIndexDat(i1,q1))
    {
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(PrimeType pt, const IQIndex& other, int inc = 1) 
	: Index(other), _dir(other._dir), pd(other.pd)  { doprime(pt,inc); }

    explicit IQIndex(istream& s) { read(s); }

    void write(ostream& s) const
    {
        Index::write(s);
        s.write((char*)&_dir,sizeof(_dir));
        pd->write(s);
    }

    void read(istream& s)
    {
        Index::read(s);
        s.read((char*)&_dir,sizeof(_dir));
        pd = new IQIndexDat(s);
    }

    IQIndex(Imaker im)
    : Index(im), _dir(In)
	{
        if(im == makeNull)
        { pd = &IQIndexDatNull; }
        else if(im == makeReIm)
        { pd = &IQIndReDat; }
        else if(im == makeReImP)
        { pd = &IQIndReDatP; }
        else if(im == makeReImPP)
        { pd = &IQIndReDatPP;}
        else if(im == makeEmptyV)
        { pd = &IQIndEmptyVDat; }
        else Error("IQIndex: Unrecognized Imaker type.");
	}

    IQIndexVal operator()(int n) const;

    //------------------------------------------
    //IQIndex: methods for querying m's

    int biggestm() const
	{
        int mm = 0;
        foreach(const inqn& x, pd->iq_) mm = max(mm,x.index.m());
        return mm;
	}
    string showm() const
	{
        string res = " ";
        std::ostringstream oh; 
        foreach(const inqn& x, pd->iq_)
        {
            QN q = x.qn;
            oh << format("[%d,%d,%s]:%d ") % q.sz() % q.Nf() % (q.fp()==1?"+":"-") % x.index.m(); 
        }
        return oh.str();
	}

    //------------------------------------------
    //IQIndex: quantum number methods

    friend inline IQIndex negate(const IQIndex& I) // Quantum numbers negated
    { 
        vector<inqn> iq(I.pd->iq_.size());
        for(size_t j = 0; j < iq.size(); ++j)
        { iq[j] = inqn(I.pd->iq_[j].index,-I.pd->iq_[j].qn); } 
        return IQIndex(I,iq); 
    }
     
    QN qn(const Index& i) const
	{ 
        foreach(const inqn& x, pd->iq_) if(x.index == i) return x.qn;
        cerr << *this << "\n";
        cerr << "i = " << i << "\n";
        Error("IQIndex::qn(Index): IQIndex does not contain given index.");
        return QN();
	}

    const Arrow dir() const { return _dir; }
    void conj() { _dir = _dir*Switch; }

    //------------------------------------------
    //IQIndex: index container methods


    const Index& findbyqn(QN q) const
	{ 
        for(size_t i = 0; i < pd->iq_.size(); ++i)
            if(pd->iq_[i].qn == q) return pd->iq_[i].index;
        Error("IQIndex::findbyqn: no Index had a matching QN.");
        return pd->iq_[0].index;
	}

    bool hasindex(const Index& i) const
	{ 
        foreach(const inqn& x, pd->iq_) 
            { if(x.index == i) return true; }
        return false;
	}
    bool hasindex_noprime(const Index& i) const
	{ 
        foreach(const inqn& x, pd->iq_) 
            { if(x.index.noprime_equals(i)) return true; }
        return false;
	}

    //------------------------------------------
    //IQIndex: prime methods

    void doprime(PrimeType pt, int inc = 1)
	{
        solo();
        Index::doprime(pt,inc);
        DoPrimer dp(pt,inc);
        for_each(pd->iq_.begin(),pd->iq_.end(),dp);
	}
    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
	{
        solo();
        Index::mapprime(plevold,plevnew,pt);
        for_each(pd->iq_.begin(),pd->iq_.end(),MapPrimer(plevold,plevnew,pt));
	}
    void noprime(PrimeType pt = primeBoth)
	{
        solo();
        Index::noprime(pt);
        for(size_t j = 0; j < pd->iq_.size(); ++j)
        { pd->iq_[j].index.noprime(pt); }
	}
    IQIndex primed(int inc = 1) const
	{
        return IQIndex(primeBoth,*this,inc);
	}

    inline friend ostream& operator<<(ostream &o, const IQIndex &I)
    {
        o << "IQIndex: " << (const Index&) I << " <" << I.dir() << ">" << endl;
        for(int j = 1; j <= I.nindex(); ++j) o << " " << I.index(j) SP I.qn(j) << "\n";
        return o;
    }

    void print(string name = "") const
    { cerr << "\n" << name << " =\n" << *this << "\n"; }

};

extern IQIndex IQIndNull, IQIndReIm, IQIndReImP, IQIndReImPP, IQEmptyV;
enum IQmaker {makeSing};

struct IQIndexVal
{
    IQIndex iqind; 
    int i;
    IQIndexVal() : iqind(IQIndNull),i(0) { }
    IQIndexVal(const IQIndex& iqindex, int i_) : iqind(iqindex),i(i_) 
    { 
        assert(i <= iqind.m());
        if(iqindex.type() != Site) Error("IQIndexVals only defined for type Site");
    }

    Index index() const { return iqind.index(i); }
    QN qn() const { return iqind.qn(i); }
    IQIndexVal primed() const { return IQIndexVal(iqind.primed(),i); }
    void conj() { iqind.conj(); }

    operator IndexVal() const { return IndexVal(iqind,i); }
    ITensor operator*(const IndexVal& iv) const { IndexVal iv_this = *this; return (iv_this * iv); }
    ITensor operator*(Real fac) const { IndexVal iv_this = *this; return iv_this * fac; }

    void print(string name = "") const
    { cerr << "\n" << name << " =\n" << *this << "\n"; }
    inline friend ostream& operator<<(ostream& s, const IQIndexVal& iv)
    { return s << "IQIndexVal: i = " << iv.i << ", iqind = " << iv.iqind << "\n"; }
};
extern IQIndexVal IQIVNull;

inline IQIndexVal IQIndex::operator()(int n) const { return IQIndexVal(*this,n); }

#ifdef THIS_IS_MAIN
IQIndexDat IQIndexDatNull(makeNull);
IQIndexDat IQIndReDat(makeReIm);
IQIndexDat IQIndReDatP(makeReImP);
IQIndexDat IQIndReDatPP(makeReImPP);
IQIndexDat IQIndEmptyVDat(makeEmptyV);

IQIndex IQIndNull(makeNull);
IQIndex IQIndReIm(makeReIm);
IQIndex IQIndReImP(makeReImP);
IQIndex IQIndReImPP(makeReImPP);
IQIndex IQEmptyV(makeEmptyV);

IQIndexVal IQIVNull(IQIndNull,0);
#endif //THIS_IS_MAIN

#endif
