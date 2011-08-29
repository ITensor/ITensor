#ifndef __IQ_H
#define __IQ_H
#include "tensor.h"
#include <set>
using std::multimap;
using std::set;

class IQCombiner;
class Condenser;

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
    int& sz() { return _sz; }
    int& Nf() { return _Nf; }
    int& Nfp() { assert(_Nfp == 0 || _Nfp == 1); return _Nfp; }

    void write(ostream& s) const { s.write((char*)this,sizeof(this)); }
    void read(istream& s) { s.read((char*)this,sizeof(this)); }

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

struct IQIndexVal;

class IQIndex : public Index
{
private:
    Arrow _dir; //Arrow direction. -1 for in, +1 for out. 
    vector<inqn> iq_;
public:
    const vector<inqn>& iq() const { return iq_; }
    const Index& index(int i) const { return GET(iq_,i-1).index; }
    const QN& qn(int i) const { return GET(iq_,i-1).qn; }

    //------------------------------------------
    //IQIndex: Constructors

    IQIndex() : _dir(Out) {}

    explicit IQIndex(const Index& other, Arrow dir = Out) : Index(other), _dir(dir) {}

    explicit IQIndex(const string& name,IndexType it = Link, Arrow dir = Out) : Index(name,1,it), _dir(dir) {}

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            Arrow dir = Out) 
    : Index(name,i1.m(),i1.type()), _dir(dir)
    {
        iq_.push_back(inqn(i1,q1));
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            Arrow dir = Out) 
    : Index(name,i1.m()+i2.m(),i1.type()), _dir(dir)
    {
        iq_.push_back(inqn(i1,q1));
        iq_.push_back(inqn(i2,q2));
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            Arrow dir = Out) 
    : Index(name,i1.m()+i2.m()+i3.m(),i1.type()), _dir(dir)
    {
        iq_.push_back(inqn(i1,q1));
        iq_.push_back(inqn(i2,q2));
        iq_.push_back(inqn(i3,q3));
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, 
            const Index& i1, const QN& q1, 
            const Index& i2, const QN& q2,
            const Index& i3, const QN& q3,
            const Index& i4, const QN& q4,
            Arrow dir = Out) 
    : Index(name,i1.m()+i2.m()+i3.m()+i4.m(),i1.type()), _dir(dir)
    {
        iq_.push_back(inqn(i1,q1));
        iq_.push_back(inqn(i2,q2));
        iq_.push_back(inqn(i3,q3));
        iq_.push_back(inqn(i4,q4));
        setPrimeLevel(i1.primelevel);
    }

    IQIndex(const string& name, const vector<inqn>& ind_qn, Arrow dir = Out) 
    : Index(name,0,ind_qn.back().index.type()), _dir(dir), iq_(ind_qn)
    { 
        int* pm = const_cast<int*>(&(p->m_));
        foreach(const inqn& x, iq_) *pm += x.index.m();
        setPrimeLevel(ind_qn.back().index.primelevel);
    }

    IQIndex(const IQIndex& other, const vector<inqn>& ind_qn)
    : Index(other.name(),0,other.type()), _dir(other._dir), iq_(ind_qn)
    { 
        int* pm = const_cast<int*>(&(p->m_));
        foreach(const inqn& x, iq_) *pm += x.index.m();
        setPrimeLevel(ind_qn.back().index.primelevel);
    }

    IQIndex(PrimeType pt, const IQIndex& other, int inc = 1) 
	: Index(other), _dir(other._dir), iq_(other.iq_)  { doprime(pt,inc); }

    IQIndex(Imaker im) : _dir(In)
	{
        Index i;
        if(im == makeNull)        { i = IndNull; }
        else if(im == makeReIm)   { i = IndReIm; }
        else if(im == makeReImP)  { i = IndReImP; }
        else if(im == makeReImPP) { i = IndReImPP; }
        else if(im == makeEmptyV) { i = IndEmptyV; }
        Index::operator=(i); 
        iq_.push_back(inqn(i,QN()));
	}

    IQIndex(istream& s) { read(s); }

    IQIndexVal operator()(int n) const;

    void write(ostream& s) const
    {
        Index::write(s);
        s.write((char*)&_dir,sizeof(_dir));
        unsigned int size = iq_.size();
        s.write((char*)&size,sizeof(size));
        foreach(const inqn& x,iq_) x.write(s);
    }

    void read(istream& s)
    {
        Index::read(s);
        s.read((char*)&_dir,sizeof(_dir));
        unsigned int size; s.read((char*)&size,sizeof(size));
        iq_.resize(size);
        foreach(inqn& x,iq_) x.read(s);
    }

    //------------------------------------------
    //IQIndex: methods for querying m's

    int biggestm() const
	{
        int mm = 0;
        foreach(const inqn& x, iq_) mm = max(mm,x.index.m());
        return mm;
	}
    string showm() const
	{
        string res = " ";
        ostringstream oh; 
        foreach(const inqn& x, iq_)
        {
            QN q = x.qn;
            oh << format("[%d,%d,%s]:%d ") % q.sz() % q.Nf() % (q.fp()==1?"+":"-") % x.index.m(); 
        }
        return oh.str();
	}

    //------------------------------------------
    //IQIndex: quantum number methods

    void negate() // negate Quantum numbers
	{ foreach(inqn& x, iq_) x.qn = -x.qn; }

    friend inline IQIndex negate(IQIndex I) // Quantum numbers negated
    { foreach(inqn& x, I.iq_) { x.qn = -x.qn; } return I; }
     
    QN qn(const Index& i) const
	{ 
        foreach(const inqn& x, iq_) if(x.index == i) return x.qn;
        cerr << *this << "\n";
        cerr << "i = " << i << "\n";
        Error("IQIndex::qn(Index): IQIndex does not contain given index.");
        return QN();
	}

    void set_qn(const Index& i, QN q)
	{
        foreach(inqn& x, iq_) if(x.index == i) { x.qn = q; return; }
        cerr << *this << "\n";
        cerr << "i = " << i << "\n";
        Error("IQIndex::qn(Index): IQIndex does not contain given index.");
	}
    void set_qn(int i, QN q)
	{ GET(iq_,i-1).qn = q; }

    const Arrow dir() const { return _dir; }
    Arrow& dir() { return _dir; }
    void conj() { _dir = _dir*Switch; }

    //------------------------------------------
    //IQIndex: index container methods

    int nindex() const { return (int) iq_.size(); }


    const Index& findbyqn(QN q) const
	{ 
        for(int i = 0; i < (int)iq_.size(); ++i)
            if(iq_[i].qn == q) return iq_[i].index;
        Error("IQIndex::findbyqn: no Index had a matching QN.");
        return iq_[0].index;
	}

    bool hasindex(const Index& i) const
	{ 
        foreach(const inqn& x, iq_) if(x.index == i) return true;
        return false;
	}
    bool hasindex_noprime(const Index& i) const
	{ 
        foreach(const inqn& x, iq_) if(x.index.noprime_equals(i)) return true;
        return false;
	}

    //------------------------------------------
    //IQIndex: prime methods

    void doprime(PrimeType pt, int inc = 1)
	{
        Index::doprime(pt,inc);
        DoPrimer dp(pt,inc);
        for_each(iq_.begin(),iq_.end(),dp);
	}
    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
	{
        Index::mapprime(plevold,plevnew,pt);
        for_each(iq_.begin(),iq_.end(),MapPrimer(plevold,plevnew,pt));
	}
    void noprime(PrimeType pt = primeBoth)
	{
        Index::noprime(pt);
        for(int j = 0; j < (int)iq_.size(); ++j)
            iq_[j].index.noprime(pt);
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
    ITensor operator*(Real fac) const { IndexVal iv_this = *this; return iv_this.operator*(fac); }
    friend inline ITensor operator*(Real fac, const IndexVal& iv) { return iv.operator*(fac); }

    void print(string name = "") const
    { cerr << "\n" << name << " =\n" << *this << "\n"; }
    inline friend ostream& operator<<(ostream& s, const IQIndexVal& iv)
    { return s << "IQIndexVal: i = " << iv.i << ", iqind = " << iv.iqind << "\n"; }
};
extern IQIndexVal IQIVNull;


class IQTensor; extern IQTensor IQTSing, IQComplex_1, IQComplex_i;

class IQTensor
{
public:
    typedef IQIndex IndexT;
    typedef IQIndexVal IndexValT;
    typedef list<ITensor>::iterator iten_it;
    typedef list<ITensor>::const_iterator const_iten_it;
    typedef vector<IQIndex>::iterator iqind_it;
    typedef vector<IQIndex>::const_iterator const_iqind_it;
private:
    bool own_rmap;
    mutable list<ITensor> itensor; // This is mutable to allow reordering
    vector<IQIndex> iqindex_;
    IQIndex viqindex; //virtual IQIndex
    map<ApproxReal,iten_it> rmap;
public:

    inline const IQIndex& iqindex(int j) const { return iqindex_[j]; }
    
    //----------------------------------------------------
    //IQTensor: iterators over 'itensor'
    iten_it       iten_begin() { return itensor.begin(); }
    iten_it       iten_end() { return itensor.end(); }
    const_iten_it const_iten_begin() const { return itensor.begin(); }
    const_iten_it const_iten_end() const { return itensor.end(); }
    pair<iten_it,iten_it> itensors() { return make_pair(iten_begin(),iten_end()); }
    pair<const_iten_it,const_iten_it> itensors() const { return make_pair(const_iten_begin(),const_iten_end()); }
    pair<iqind_it,iqind_it> iqinds() { return make_pair(iqindex_.begin(),iqindex_.end()); }
    pair<const_iqind_it,const_iqind_it> iqinds() const { return make_pair(iqindex_.begin(),iqindex_.end()); }

    //----------------------------------------------------
    //IQTensor: Constructors

    IQTensor() : own_rmap(false), viqindex(IQEmptyV) {}


    explicit IQTensor(const IQIndex& i1) :  own_rmap(false), viqindex(IQEmptyV)
    { iqindex_.push_back(i1); }
    IQTensor(const IQIndex& i1,const IQIndex& i2) : own_rmap(false), viqindex(IQEmptyV)
    { iqindex_.push_back(i1); iqindex_.push_back(i2); }
    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3) : own_rmap(false), viqindex(IQEmptyV)
    { iqindex_.push_back(i1); iqindex_.push_back(i2); iqindex_.push_back(i3); }
    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,const IQIndex& i4) : own_rmap(false), viqindex(IQEmptyV)
    { iqindex_.push_back(i1); iqindex_.push_back(i2); iqindex_.push_back(i3); iqindex_.push_back(i4); }

    explicit IQTensor(const vector<IQIndex>& iqinds_) : own_rmap(false), viqindex(IQEmptyV)
    { iqindex_ = iqinds_; }

    IQTensor(ITmaker itm) : own_rmap(false), viqindex(IQEmptyV)
    {
        iqindex_.push_back(IQIndReIm);
        if(itm == makeComplex_1) 
            operator+=(Complex_1);
        if(itm == makeComplex_i) 
            operator+=(Complex_i);
    }

    IQTensor(IQmaker i) : own_rmap(false), viqindex(IQEmptyV)
    {
        Index s("sing");
        IQIndex single("single",s,QN());
        iqindex_.push_back(single);
        ITensor st(s);
        st.ncdat() = 1.0;
        operator+=(st);
    }

    IQTensor(PrimeType pt,const IQTensor& other) 
    : own_rmap(false), itensor(other.itensor), iqindex_(other.iqindex_), viqindex(other.viqindex)
    { doprime(pt); }

    IQTensor(const IQTensor& other)
    : own_rmap(false), itensor(other.itensor), iqindex_(other.iqindex_), viqindex(other.viqindex)  
    { }

    IQTensor& operator=(const IQTensor& other)
    {
        if(this == &other) return *this;
        iqindex_ = other.iqindex_;
        viqindex = other.viqindex;
        itensor = other.itensor;
        rmap.clear();
        own_rmap = false;
        return *this;
    }

    //----------------------------------------------------
    //IQTensor operators

    IQTensor operator*(IQTensor other) const { other *= *this; return other; }
    IQTensor& operator*=(const IQTensor& other);

    IQTensor& operator+=(const IQTensor& o);
    IQTensor operator+(const IQTensor& o) const { IQTensor res(*this); res += o; return res; }
    IQTensor operator-(const IQTensor& o) const 
    { IQTensor res(o); res *= -1.0; res += *this; return res; }

    IQTensor& operator*=(Real fac) 
    { 
        if(fac == 0.0) { itensor.clear(); own_rmap = false;  return *this; }
        foreach(ITensor& t, itensor) { t *= fac; }
        return *this; 
    }
    IQTensor operator*(Real fac) { IQTensor res(*this); res *= fac; return res; }

    /*
    operator ITensor() const
    {
        vector<Index> indices;
        foreach(const IQIndex& I, iqindex_)
        {
            if(I.type() != Site) 
            { Error("IQTensor to ITensor conversion requires all IQIndex's of type Site."); }
            indices.push_back(I);
        }
        const bool do_allocate = true;
        ITensor res(indices,do_allocate);
        match_order();
        return res;
    }
    */

    void insert(const ITensor& t) 
    { 
        if(!own_rmap)
        {
            rmap.clear();
            iten_it b = iten_begin();
            iten_it e = iten_end();
            for(iten_it jj = b; jj != e; ++jj)
                rmap.insert(make_pair(ApproxReal(jj->unique_Real()),jj));
            own_rmap = true;
        }
        ApproxReal r(t.unique_Real());
        if(rmap.count(r) == 0)
        {
            itensor.push_front(t);
            rmap.insert(pair<ApproxReal,iten_it>(r,itensor.begin()));
        }
        else
        {
            cerr << "t is " << t << endl;
            Error("Can't insert ITensor with identical structure twice.");
        }
    } //end IQTensor::insert(const ITensor& t)

    IQTensor& operator+=(const ITensor& t) 
    { 
        if(!own_rmap)
        {
            rmap.clear();
            iten_it b = iten_begin();
            iten_it e = iten_end();
            for(iten_it jj = b; jj != e; ++jj)
            { rmap.insert(make_pair(ApproxReal(jj->unique_Real()),jj)); }
            own_rmap = true;
        }
        ApproxReal r(t.unique_Real());
        if(rmap.count(r) == 0)
        {
            itensor.push_front(t);
            rmap.insert(pair<ApproxReal,iten_it>(r,itensor.begin()));
        }
        else *(rmap[r]) += t;
        return *this;
    } //end IQTensor::operator+=(ITensor)

    Real& operator()(const IQIndexVal& iv1, const IQIndexVal& iv2 = IQIVNull, const IQIndexVal& iv3 = IQIVNull,
                     const IQIndexVal& iv4 = IQIVNull, const IQIndexVal& iv5 = IQIVNull, const IQIndexVal& iv6 = IQIVNull,
                     const IQIndexVal& iv7 = IQIVNull, const IQIndexVal& iv8 = IQIVNull)
	{
        array<IQIndexVal,NMAX+1> iv = {{ IQIVNull, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};
        Real ur = 0; int nn = 0; 
        while(iv[nn+1].iqind != IQIVNull.iqind) 
        { ++nn; ur += GET(iv,nn).index().unique_Real(); }
        ApproxReal r(ur);

        if(rmap.count(r) == 0)
        {
            vector<Index> indices; indices.reserve(nn);
            foreach(const IQIndex& I, iqindex_)
            {
                if(I.type() == Site) continue;
                if(I.m() == 1 && I.nindex() == 1) indices.push_back(I.index(1));
                else Error("IQTensor::operator() not permitted for IQTensor with non-trivial Link indices.");
            }
            for(int j = 1; j <= nn; ++j) 
            {
                if(!hasindex(iv[j].iqind)) Error("IQTensor::operator(): IQIndex not found.");
                indices.push_back(iv[j].index());
            }
            ITensor t(indices,true);
            itensor.push_front(t);
            rmap.insert(pair<ApproxReal,iten_it>(r,itensor.begin()));
        }
        return (rmap[r])->operator()(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8);
    }

    //----------------------------------------------------
    //IQTensor quantum number methods
    void set_qn(Index i, QN q)
    {
        int iqq = find_iqind(i);
        if(iqq == -1)
            Error("set_qn: cant find index");
        iqindex_[iqq].set_qn(i,q);
    } //end IQTensor::set_qn

    QN net_QN() const // only works on tensors without the Link indices put in
    {
        assert(!hastype(Link));
        QN res;
        for(const_iten_it t = const_iten_begin(); t != const_iten_end(); t++)
        {
            QN q;
            for(int j = 1; j <= t->r(); j++)
            { q += qn(t->index(j))*dir(t->index(j)); }
            
            if(t == const_iten_begin()) res = q;
            else if(q != res)
            {
                printdat = true; cerr << "*this = " << *this << endl; printdat = false;
                cerr << "res = " << res << endl;
                cerr << "q = " << q << endl;
                error("net_QN() should only be used on an IQTensor with well defined QN.");
            }
        }
        return res;
    } //end IQTensor::net_QN

    QN qn(const Index& in) const
    {
        int iqq = find_iqind(in);
        if(iqq == -1)
            Error("qn: cant find index");
        return iqindex_[iqq].qn(in);
    } //end IQTensor::qn

    void negate_qn() { foreach(IQIndex& I, iqindex_) I.negate(); }

    Arrow dir(const Index& in) const
    {
        int iqq = find_iqind(in);
        if(iqq == -1) 
        {
            this->print("this IQTensor");
            in.print("in"); 
            Error("IQTensor::dir(Index&): cant find Index in IQIndices");
        }
        return iqindex_[iqq].dir();
    } //end IQTensor::dir

    //Checks if the divergence of this IQTensor is zero
    bool checkDivZero() const
    {
        foreach(const ITensor& it, itensor)
        {
            QN qtot;
            for(int j = 1; j <= it.r(); ++j) 
            { qtot += qn(it.index(j))*dir(it.index(j)); }
            qtot += viqindex.qn(1)*viqindex.dir();
            if(qtot != QN()) 
            {
                cerr << "checkDivZero: IQTensor failed to have zero divergence." << endl;
                cerr << "qtot = " << qtot << "\n";
                printdat = true; cerr << "Offending ITensor = " << it << "\n"; printdat = false;
                cerr << "*this = "; print();
                return false;
            }
        }
        return true;
    }

    //----------------------------------------------------
    //IQTensor: prime methods

    void ind_inc_prime(const IQIndex& i,int inc)
    {
        own_rmap = false;
        bool gotit = false;
        foreach(IQIndex& jj, iqindex_)
        if(jj.noprime_equals(i))
        {
            gotit = true;
            int p = jj.primelevel;
            jj.mapprime(p,p+inc);
        }

        if(!gotit)
        {
            cerr << "IQIndex was " << i << "\n";
            Error("ind_inc_prime: couldn't find IQIndex");
        }

        for(iten_it jj = iten_begin(); jj != iten_end(); ++jj)
        for(int ii = 1; ii <= jj->r(); ++ii)
        if(i.hasindex_noprime(jj->index(ii)))
        {
            int p = jj->index(ii).primelevel;
            jj->mapprimeind(jj->index(ii),p,p+inc);
        }
    } //end IQTensor::ind_inc_prime

    void noprime(PrimeType pt = primeBoth)
    {
        own_rmap = false;
        foreach(IQIndex& J, iqindex_)
	    J.noprime(pt);

        for(iten_it jj = iten_begin(); jj != iten_end(); ++jj)
	    foreach(ITensor& t, itensors()) t.noprime(pt);
    } //end IQTensor::noprime
    friend inline IQTensor deprimed(IQTensor A)
    { A.noprime(); return A; }

    void noprimelink()
    {
        own_rmap = false;
        foreach(IQIndex& J, iqindex_)
	    if(J.type() == Link) J.noprime();

        for(iten_it jj = iten_begin(); jj != iten_end(); ++jj)
	    foreach(ITensor& t, itensors()) t.noprime(primeLink);
    } //end IQTensor::noprimelink

    void doprime(PrimeType pt)
    {
        own_rmap = false;
        DoPrimer prim(pt);
        for_each(iqindex_.begin(), iqindex_.end(),prim);
        for(iten_it jj = iten_begin(); jj != iten_end(); ++jj)
        jj->doprime(pt);
    } //end IQTensor::doprime

    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth) // no need to keep prime level small
    {
        own_rmap = false;
        MapPrimer prim(plevold,plevnew,pt);
        for_each(iqindex_.begin(), iqindex_.end(),prim);
        for(iten_it jj = iten_begin(); jj != iten_end(); ++jj)
        jj->mapprime(plevold,plevnew,pt);
    } //end IQTensor::mapprime

    void primeind(const IQIndex& I)
    {
        foreach(IQIndex& J, iqindex_)
        { if(J == I) J = J.primed(); }
        if(viqindex == I) viqindex = viqindex.primed();

        foreach(ITensor& t, itensor)
        foreach(const inqn& x, I.iq())
        { if(t.hasindex(x.index)) t.primeind(x.index); }
    }
    friend inline IQTensor primeind(IQTensor A, const IQIndex& I)
    { A.primeind(I); return A; }

    friend inline IQTensor primed(IQTensor A)
    { A.doprime(primeBoth); return A; }

    void primesite() { doprime(primeSite); }
    friend inline IQTensor primesite(IQTensor A)
    { A.doprime(primeSite); return A; }

    void primelink() { doprime(primeLink); }
    friend inline IQTensor primelink(const IQTensor& A)
    { IQTensor res(A); res.doprime(primeLink); return res; }


    //----------------------------------------------------
    //IQTensor index methods

    int find_iqind(const Index& I) const
    {
        for(unsigned int j = 0; j < iqindex_.size(); ++j)
        if(iqindex_[j].hasindex(I)) { return j; }
        return -1;
    }

    //Return true if one of the ITensors uses this Index
    bool uses_ind(const Index& i) const
    {
        foreach(const ITensor& it, itensor)
        if(it.hasindex(i)) { return true; }
        return false;
    }

    int findindex(const IQIndex& i) const
    {
        vector<IQIndex>::const_iterator f = find(iqindex_.begin(),iqindex_.end(),i);
        if(f == iqindex_.end()) return -1;
        else return (f - iqindex_.begin());
    }

    bool hastype(IndexType t) const
    {
        if(t == Virtual) return (viqindex != IQEmptyV);
        foreach(const IQIndex& I, iqindex_)
        if(I.type() == t) { return true; }
        return false;
    }

    const IQIndex& findtype(IndexType t) const
    {
        if(t == Virtual)
        {
            if(viqindex != IQEmptyV) return viqindex;
            else Error("IQTensor::findtype: couldn't find an IQIndex of type Virtual.");
        }
        foreach(const IQIndex& I, iqindex_)
        if(I.type() == t) { return I; }
        Error("IQTensor::findtype: couldn't find type");
        return viqindex;
    }

    const IQIndex& finddir(Arrow dir) const
    {
        foreach(const IQIndex& I, iqindex_)
        if(I.dir() == dir) { return I; }
        Error("IQTensor::finddir: couldn't find dir");
        return viqindex;
    }

    bool hasindex(const IQIndex& i) const { return findindex(i) != -1; }
    bool has_virtual() const { return viqindex != IQEmptyV; }
    QN virtualQN() const { return viqindex.qn(1); }
    const IQIndex& virtual_ind() const { return viqindex; }

    void addindex1(IQIndex I)
    {
        if(I.m() != 1) Error("IQTensor::operator*=(IQIndex): IQIndex must have m == 1.");    
        if(I.type() == Virtual) 
        {
            if(viqindex != IQEmptyV) //Add quantum numbers (with arrows)
            {
                QN newq = I.dir()*(viqindex.dir()*viqindex.qn(1)+I.dir()*I.qn(1));
                if(newq.Nf() < 0) //prefer to have Nf >= 0
                {
                    newq *= -1;
                    I.conj();
                }
                I.set_qn(1,newq);
            }
            viqindex = I;
            return;
        }
        foreach(ITensor& t, itensor) t.addindex1(I.index(1));
        iqindex_.push_back(I);
    }

    //----------------------------------------------------
    //IQTensor miscellaneous methods

    Real unique_Real() const
    {
        Real ur = 0.0;
        foreach(const IQIndex& I, iqindex_) ur += I.unique_Real();
        return ur;
    }
    inline int iten_size() const { return itensor.size(); }
    inline bool is_null() const { return iten_size() == 0; }
    inline bool is_not_null() const { return iten_size() != 0; }
    inline int num_index() const { return iqindex_.size(); }
    int num_index(IndexType t) const 
    { 
        int count = 0;
        foreach(const IQIndex& I, iqindex_)
        if(I.type() == t) ++count;
        return count;
    }

    Real norm() const
    {
        Real res;
        foreach(const ITensor& t, itensor) res += sqrt(t.norm());
        return sqrt(res);
    }

    int vec_size() const
    {
        int s = 0;
        for(const_iten_it jj = const_iten_begin(); jj != const_iten_end(); ++jj)
            s += jj->Length();
        return s;
    }
    void AssignToVec(VectorRef v) const
    {
        if(vec_size() != v.Length())
            Error("Mismatched sizes in IQTensor::AssignToVec(VectorRef v).");
        int off = 1;
        for(const_iten_it jj = const_iten_begin(); jj != const_iten_end(); ++jj)
            {
            int d = jj->Length();
            jj->AssignToVec(v.SubVector(off,off+d-1));
            off += d;
            }
    }
    void AssignFromVec(VectorRef v)
    {
        if(vec_size() != v.Length())
            Error("bad size");
        int off = 1;
        for(iten_it jj = iten_begin(); jj != iten_end(); ++jj)
            {
            int d = jj->dat().Length();
            jj->AssignFromVec(v.SubVector(off,off+d-1));
            off += d;
            }
    }
    void GetSingComplex(Real& re, Real& im) const;

    
    void Randomize() { foreach(ITensor& t, itensor) t.Randomize(); }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    void Assign(const IQTensor& other);

    void SplitReIm(IQTensor& re, IQTensor& im) const;

    friend IQTensor conj(const IQTensor& A)
    {
        IQTensor r,i;
        if(!A.hasindex(IQIndReIm))
        {
            IQTensor r = A;
            foreach(IQIndex& I,r.iqindex_) I.conj();
            return r;
        }
        else
        {
            A.SplitReIm(r,i);
            foreach(IQIndex& I,r.iqindex_) I.conj();
            foreach(IQIndex& I,i.iqindex_) I.conj();
            i *= -1.0;
            return r * IQComplex_1 + IQComplex_i * i;
        }
    }

private:
    void match_order() const
    {
        const int s = iqindex_.size();
        foreach(ITensor& t, itensor)
        {
            if(t.r() != s)
            {
                Print(*this); Print(t);
                Error("match_order: ds not same");
            }
            Permutation P;
            for(int j = 1; j <= t.r(); ++j)
            {
                bool gotone = false;
                for(int k = 0; k < s; ++k)
                if(iqindex_[k].hasindex(t.index(j)))
                { P.from_to(j,k+1); gotone = true; break; }
                if(!gotone) Error("match_order: !gotone");
            }
            t.Reshape(P);
        }
    }

public:
    friend IQTensor operator*(const IQTensor& t, const IQCombiner& c);
    friend IQTensor operator*(const IQTensor& t, const Condenser& c);

    inline friend ostream& operator<<(ostream & s, const IQTensor &t)
    {
        s << "\n----- IQTensor -----\nIQIndices: " << endl;
        foreach(const IQIndex& I, t.iqinds()) s << "  " << I << endl;
        if(t.has_virtual()) s << "  " << t.virtual_ind() << endl;
        s << "ITensors: \n";
        foreach(const ITensor& i, t.itensors())
        s <<"	" << i << "\n";
        s << "-------------------" << "\n\n";
        return s;
    }

}; //class IQTensor

IQIndex index_in_common(const IQTensor& A, const IQTensor& B, IndexType t);

inline Real ReSingVal(const IQTensor& x)
{
    Real re, im;
    x.GetSingComplex(re,im);
    return re;
}
inline Real Dot(const IQTensor& x, const IQTensor& y, bool doconj = true)
{
    IQTensor res(IQTSing*(doconj ? conj(x) : x)*y);
    return ReSingVal(res);
}
inline void Dot(const IQTensor& x, const IQTensor& y, Real& re, Real& im, bool doconj = true)
{
    IQTensor res(IQTSing*(doconj ? conj(x) : x)*y);
    res.GetSingComplex(re,im);
}

inline IQTensor operator*(Real fac, IQTensor T)
    { T *= fac; return T; }

template<class T, class C>
bool has_element(const T& t, const C& c)
{
    return find(c.begin(),c.end(),t) != c.end();
}


template<class T> class Printit
{
public:
    ostream& s;
    string spacer;
    Printit(ostream& _s, string _spacer) : s(_s), spacer(_spacer) {}
    void operator()(const T& t) { s << t << spacer; }
};


class SiteOp
{
    const IQIndex si;
    mutable bool made_iqt, made_t;
    mutable IQTensor iqt;
    mutable ITensor t;
    map<ApproxReal, pair<IQIndexVal,IQIndexVal> > ivmap;
    map<ApproxReal, Real> valmap;
    typedef map<ApproxReal, Real>::value_type valmap_vt;

    pair<IQIndexVal,IQIndexVal> civmap(const ApproxReal& key) const { return ivmap.find(key)->second; }

    void make_iqt() const
    {
        if(made_iqt) return;
        iqt = IQTensor(conj(si),si.primed());
        foreach(const valmap_vt& x, valmap)
        { iqt(civmap(x.first).first,civmap(x.first).second) = x.second; }
        made_iqt = true;
    }
    void make_t() const
    {
        if(made_t) return;
        t = ITensor(conj(si),si.primed());
        foreach(const valmap_vt& x, valmap)
        { t(civmap(x.first).first,civmap(x.first).second) = x.second; }
        made_t = true;
    }
public:
    SiteOp(const IQIndex& si_) : si(si_), made_iqt(false), made_t(false) { }
    operator ITensor() const { make_t(); return t; }
    operator IQTensor() const { make_iqt(); return iqt; }

    Real& operator()(const IQIndexVal& iv, const IQIndexVal& ivp)
    {
        if(iv.iqind.primelevel != 0) Error("SiteOp::operator(): first IndexVal must be unprimed.");
        if(ivp.iqind.primelevel != 1) Error("SiteOp::operator(): second IndexVal must be primed.");
        Real r = iv.iqind.unique_Real() + ivp.iqind.unique_Real() + iv.i + 1000*ivp.i;
        ivmap[ApproxReal(r)] = make_pair(iv,ivp);
        return valmap[ApproxReal(r)];
    }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }
    friend ostream& operator<<(ostream& s, const SiteOp& op) { return s << IQTensor(op); }

    template <class X>
    ITensor operator*(const X& x) const { ITensor res(*this); res *= x; return res; }
    template <class X>
    friend inline ITensor operator*(const X& x, const SiteOp& op) { ITensor res(op); return (res *= x); }

    IQTensor operator*(const IQTensor& t) const { IQTensor res(*this); res *= t; return res; }
    friend inline IQTensor operator*(const IQTensor& t, const SiteOp& op) { IQTensor res(op); return (res *= t); }

};

#ifdef THIS_IS_MAIN
IQIndex IQIndNull(makeNull);
IQIndex IQIndReIm(makeReIm);
IQIndex IQIndReImP(makeReImP);
IQIndex IQIndReImPP(makeReImPP);
IQIndex IQEmptyV(makeEmptyV);
IQTensor IQTSing(makeSing);
IQTensor IQComplex_1(makeComplex_1), IQComplex_i(makeComplex_i);
IQIndexVal IQIVNull(IQIndNull,0);
#endif
#endif
