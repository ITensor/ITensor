#ifndef __IQ_H
#define __IQ_H
#include "itensor.h"
#include <set>
#include <list>
#include <map>
using std::list;
using std::map;

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

    int nindex() const { return (int) pd->iq_.size(); }

    const Index& findbyqn(QN q) const
	{ 
        for(size_t i = 0; i < pd->iq_.size(); ++i)
            if(pd->iq_[i].qn == q) return pd->iq_[i].index;
        Error("IQIndex::findbyqn: no Index had a matching QN.");
        return pd->iq_[0].index;
	}

    bool hasindex(const Index& i) const
	{ 
        foreach(const inqn& x, pd->iq_) if(x.index == i) return true;
        return false;
	}
    bool hasindex_noprime(const Index& i) const
	{ 
        foreach(const inqn& x, pd->iq_) if(x.index.noprime_equals(i)) return true;
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

class IQTDat
{
    typedef list<ITensor>::iterator       iten_it;
    typedef list<ITensor>::const_iterator const_iten_it;
private:
    mutable unsigned int numref;
    mutable bool rmap_init;
public:
    mutable list<ITensor> itensor; // This is mutable to allow reordering
    vector<IQIndex> iqindex_;
    mutable map<ApproxReal,iten_it> rmap; //mutable so that const IQTensor methods can use rmap

    IQTDat() : numref(0), rmap_init(false) { }

    explicit IQTDat(const IQIndex& i1): numref(0), rmap_init(false), iqindex_(1)
    { 
        assert(i1.type() != Virtual); iqindex_[0] = i1; 
    }

    IQTDat(const IQIndex& i1, const IQIndex& i2): numref(0), rmap_init(false), iqindex_(2)
    { 
        assert(i1.type() != Virtual); iqindex_[0] = i1; 
        assert(i2.type() != Virtual); iqindex_[1] = i2; 
    }

    IQTDat(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3): numref(0), rmap_init(false), iqindex_(3)
    { 
        assert(i1.type() != Virtual); iqindex_[0] = i1; 
        assert(i2.type() != Virtual); iqindex_[1] = i2; 
        assert(i3.type() != Virtual); iqindex_[2] = i3; 
    }

    IQTDat(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3, const IQIndex& i4): numref(0), rmap_init(false), iqindex_(4)
    { 
        assert(i1.type() != Virtual); iqindex_[0] = i1; 
        assert(i2.type() != Virtual); iqindex_[1] = i2; 
        assert(i3.type() != Virtual); iqindex_[2] = i3; 
        assert(i4.type() != Virtual); iqindex_[3] = i4; 
    }

    explicit IQTDat(vector<IQIndex>& iqinds_) : numref(0), rmap_init(false) { iqindex_.swap(iqinds_); }

    explicit IQTDat(const IQTDat& other) 
    : numref(0), rmap_init(false), itensor(other.itensor), iqindex_(other.iqindex_)
    { }

    explicit IQTDat(istream& s) : numref(0) { read(s); }

    void read(istream& s)
    {
        rmap_init = false;
        size_t size;
        s.read((char*) &size,sizeof(size));
        itensor.resize(size);
        foreach(ITensor& t, itensor) t.read(s);

        s.read((char*) &size,sizeof(size));
        iqindex_.resize(size);
        foreach(IQIndex& I, iqindex_) I.read(s);
    }

    void write(ostream& s) const
    {
        size_t size = itensor.size();
        s.write((char*) &size,sizeof(size));
        foreach(const ITensor& t, itensor) t.write(s);

        size = iqindex_.size();
        s.write((char*) &size,sizeof(size));
        foreach(const IQIndex& I, iqindex_) I.write(s);
    }

    void init_rmap() const
    {
        if(rmap_init) return;
        for(iten_it it = itensor.begin(); it != itensor.end(); ++it)
        { rmap[ApproxReal(it->unique_Real())] = it; }
        rmap_init = true;
    }

    void uninit_rmap() const { assert(numref == 1); rmap_init = false; }

    ENABLE_INTRUSIVE_PTR(IQTDat)
private:
    ~IQTDat() { } //must be dynamically allocated
    void operator=(const IQTDat&);
};

class IQTensor; extern IQTensor IQTSing, IQComplex_1, IQComplex_i;

class IQTensor
{
public:
    typedef IQIndex IndexT;
    typedef IQIndexVal IndexValT;
    typedef IQCombiner CombinerT;
    typedef list<ITensor>::iterator iten_it;
    typedef list<ITensor>::const_iterator const_iten_it;
    typedef vector<IQIndex>::iterator iqind_it;
    typedef vector<IQIndex>::const_iterator const_iqind_it;
    static const IQIndex& ReImIndex;
private:
    intrusive_ptr<IQTDat> p;
    IQIndex viqindex; //virtual IQIndex

    void solo()
    {
        assert(p != 0);
        if(p->count() != 1)
        {
            intrusive_ptr<IQTDat> new_p(new IQTDat(*p));
            p.swap(new_p);
        }
    }
public:
    int r() const { return p->iqindex_.size(); }
    inline const IQIndex& index(int j) const { return GET(p->iqindex_,j-1); }
    inline int iten_size() const { return p->itensor.size(); }
    inline bool is_null() const { return p == 0; }
    inline bool is_not_null() const { return p != 0; }
    inline int num_index() const { return p->iqindex_.size(); }
    
    //----------------------------------------------------
    //IQTensor: iterators 
    const_iten_it const_iten_begin() const { return p->itensor.begin(); }
    const_iten_it const_iten_end() const { return p->itensor.end(); }
    pair<const_iten_it,const_iten_it> itensors() const { return make_pair(p->itensor.begin(),p->itensor.end()); }

    const_iqind_it const_iqind_begin() const { return p->iqindex_.begin(); }
    const_iqind_it const_iqind_end()   const { return p->iqindex_.end(); }
    pair<const_iqind_it,const_iqind_it> iqinds() const { return make_pair(p->iqindex_.begin(),p->iqindex_.end()); }

    //----------------------------------------------------
    //IQTensor: Constructors

    IQTensor() : p(0), viqindex(IQEmptyV) {}

    explicit IQTensor(const IQIndex& i1) 
    : p(new IQTDat(i1)), viqindex(IQEmptyV)
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2) 
    : p(new IQTDat(i1,i2)), viqindex(IQEmptyV)
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3) 
    : p(new IQTDat(i1,i2,i3)), viqindex(IQEmptyV)
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,const IQIndex& i4) 
    : p(new IQTDat(i1,i2,i3,i4)), viqindex(IQEmptyV)
    { }

    explicit IQTensor(vector<IQIndex>& iqinds_) 
    : p(new IQTDat(iqinds_)), viqindex(IQEmptyV)
    { }

    IQTensor(ITmaker itm) : p(new IQTDat(IQIndReIm)), viqindex(IQEmptyV)
    {
        if(itm == makeComplex_1)      operator+=(Complex_1);
        else if(itm == makeComplex_i) operator+=(Complex_i);
    }

    IQTensor(IQmaker i) : p(new IQTDat()), viqindex(IQEmptyV)
    {
        Index s("sing");
        IQIndex single("single",s,QN());
        p->iqindex_.push_back(single);
        ITensor st(s,1);
        operator+=(st);
    }

    IQTensor(PrimeType pt,const IQTensor& other) : p(other.p), viqindex(other.viqindex)
    { doprime(pt); }

    explicit IQTensor(istream& s) : p(0), viqindex(IQEmptyV) { read(s); }

    void read(istream& s)
    {
        bool null_;
        s.read((char*) &null_,sizeof(null_));
        if(null_) { *this = IQTensor(); return; }
        p = new IQTDat(s);
        bool has_v;
        s.read((char*) &has_v,sizeof(has_v));
        if(has_v) viqindex.read(s);
    }

    void write(ostream& s) const
    {
        bool null_ = is_null();
        s.write((char*) &null_,sizeof(null_));
        if(null_) return;
        p->write(s);
        bool has_v = has_virtual();
        s.write((char*) &has_v,sizeof(has_v));
        if(has_v) viqindex.write(s);
    }


    //----------------------------------------------------
    //IQTensor operators

    IQTensor operator*(IQTensor other) const { other *= *this; return other; }
    IQTensor& operator*=(const IQTensor& other);

    IQTensor& operator+=(const IQTensor& o);
    IQTensor operator+(const IQTensor& o) const { IQTensor res(*this); res += o; return res; }
    IQTensor operator-(const IQTensor& o) const 
    { IQTensor res(o); res *= -1; res += *this; return res; }

    IQTensor& operator*=(Real fac) 
    { 
        solo();
        if(fac == 0) { p->itensor.clear(); p->uninit_rmap(); return *this; }
        foreach(ITensor& t, p->itensor) { t *= fac; }
        return *this; 
    }
    IQTensor operator*(Real fac) { IQTensor res(*this); res *= fac; return res; }
    friend inline IQTensor operator*(Real fac, IQTensor T) { T *= fac; return T; }

    /*
    operator ITensor() const
    {
        vector<Index> indices;
        foreach(const IQIndex& I, p->iqindex_)
        {
            if(I.type() != Site) 
            { Error("IQTensor to ITensor conversion requires all IQIndex's of type Site."); }
            indices.push_back(I);
        }
        ITensor res(indices);
        match_order();
        return res;
    }
    */

    void insert(const ITensor& t) 
    { 
        solo();
        p->init_rmap();
        ApproxReal r(t.unique_Real());
        if(p->rmap.count(r) == 0)
        {
            p->itensor.push_front(t);
            p->rmap.insert(pair<ApproxReal,iten_it>(r,p->itensor.begin()));
            return;
        }
        Print(*(p->rmap[r])); Print(t);
        Error("IQTensor:: Can't insert ITensor with identical structure twice, use operator+=.");
    } //IQTensor::insert(const ITensor& t)

    IQTensor& operator+=(const ITensor& t) 
    { 
        solo();
        p->init_rmap();
        ApproxReal r(t.unique_Real());
        if(p->rmap.count(r) == 0)
        {
            p->itensor.push_front(t);
            p->rmap[r] = p->itensor.begin();
        }
        else *(p->rmap[r]) += t;
        return *this;
    } //IQTensor::operator+=(ITensor)

    Real& operator()(const IQIndexVal& iv1, const IQIndexVal& iv2 = IQIVNull, const IQIndexVal& iv3 = IQIVNull,
                     const IQIndexVal& iv4 = IQIVNull, const IQIndexVal& iv5 = IQIVNull, const IQIndexVal& iv6 = IQIVNull,
                     const IQIndexVal& iv7 = IQIVNull, const IQIndexVal& iv8 = IQIVNull)
	{
        solo(); p->init_rmap();
        array<IQIndexVal,NMAX+1> iv = {{ IQIVNull, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};
        Real ur = 0; int nn = 0; 
        while(GET(iv,nn+1).iqind != IQIVNull.iqind) 
        { ++nn; ur += GET(iv,nn).index().unique_Real(); }
        ApproxReal r(ur);

        if(p->rmap.count(r) == 0)
        {
            vector<Index> indices; indices.reserve(nn);
            foreach(const IQIndex& I, p->iqindex_)
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
            p->itensor.push_front(t);
            p->rmap[r] = p->itensor.begin();
        }
        return (p->rmap[r])->operator()(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8);
    }

    //----------------------------------------------------
    //IQTensor quantum number methods
    /*
    void set_qn(Index i, QN q)
    {
        int iqq = find_iqind(i)-1;
        if(iqq == -1)
            Error("set_qn: cant find index");
        p->iqindex_[iqq].set_qn(i,q);
    } //end IQTensor::set_qn
    */

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
        int iqq = find_iqind(in)-1;
        if(iqq == -1) Error("qn: cant find index");
        return p->iqindex_[iqq].qn(in);
    } //end IQTensor::qn

    //void negate_qn() { foreach(IQIndex& I, p->iqindex_) I.negate(); }

    Arrow dir(const Index& in) const
    {
        int iqq = find_iqind(in)-1;
        if(iqq == -1) 
        {
            this->print("this IQTensor");
            in.print("in"); 
            Error("IQTensor::dir(Index&): cant find Index in IQIndices");
        }
        return p->iqindex_[iqq].dir();
    } //end IQTensor::dir


    //----------------------------------------------------
    //IQTensor: prime methods

    void ind_inc_prime(const IQIndex& i,int inc)
    {
        if(viqindex == i) { int p = viqindex.primelevel; viqindex.mapprime(p,p+inc); return; }

        solo();
        p->uninit_rmap();
        bool gotit = false;
        foreach(IQIndex& jj, p->iqindex_)
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

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
        for(int ii = 1; ii <= jj->r(); ++ii)
        if(i.hasindex_noprime(jj->index(ii)))
        {
            int p = jj->index(ii).primelevel;
            jj->mapprimeind(jj->index(ii),p,p+inc);
        }
    } //end IQTensor::ind_inc_prime

    void noprime(PrimeType pt = primeBoth)
    {
        solo();
        p->uninit_rmap();
        foreach(IQIndex& J, p->iqindex_)
	    J.noprime(pt);

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    foreach(ITensor& t, p->itensor) t.noprime(pt);
    } //end IQTensor::noprime
    friend inline IQTensor deprimed(IQTensor A)
    { A.noprime(); return A; }

    void noprimelink()
    {
        solo();
        p->uninit_rmap();
        foreach(IQIndex& J, p->iqindex_)
	    if(J.type() == Link) J.noprime();

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    foreach(ITensor& t, p->itensor) t.noprime(primeLink);
    } //end IQTensor::noprimelink

    void doprime(PrimeType pt)
    {
        solo();
        p->uninit_rmap();
        DoPrimer prim(pt);
        for_each(p->iqindex_.begin(), p->iqindex_.end(),prim);
        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
        { jj->doprime(pt); }
        viqindex.doprime(pt);
    } //end IQTensor::doprime

    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth) // no need to keep prime level small
    {
        solo();
        p->uninit_rmap();
        MapPrimer prim(plevold,plevnew,pt);
        for_each(p->iqindex_.begin(), p->iqindex_.end(),prim);
        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
        jj->mapprime(plevold,plevnew,pt);
        viqindex.mapprime(plevold,plevnew,pt);
    } //end IQTensor::mapprime

    void primeind(const IQIndex& I)
    {
        if(viqindex == I) { viqindex = viqindex.primed(); return; }

        solo();
        p->uninit_rmap();
        foreach(IQIndex& J, p->iqindex_)
        { if(J == I) J = J.primed(); }

        foreach(ITensor& t, p->itensor)
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
        for(size_t j = 0; j < p->iqindex_.size(); ++j)
        if(p->iqindex_[j].hasindex(I)) { return j+1; }
        return 0;
    }

    //Return true if one of the ITensors uses this Index
    bool uses_ind(const Index& i) const
    {
        foreach(const ITensor& it, p->itensor)
        if(it.hasindex(i)) { return true; }
        return false;
    }

    int findindex(const IQIndex& i) const
    {
        vector<IQIndex>::const_iterator f = find(p->iqindex_.begin(),p->iqindex_.end(),i);
        if(f == p->iqindex_.end()) return 0;
        else return (f - p->iqindex_.begin())+1;
    }

    bool hastype(IndexType t) const
    {
        if(t == Virtual) return (viqindex != IQEmptyV);
        foreach(const IQIndex& I, p->iqindex_)
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
        foreach(const IQIndex& I, p->iqindex_)
        if(I.type() == t) { return I; }
        Error("IQTensor::findtype: couldn't find type");
        return viqindex;
    }

    const IQIndex& finddir(Arrow dir) const
    {
        foreach(const IQIndex& I, p->iqindex_)
        if(I.dir() == dir) { return I; }
        Error("IQTensor::finddir: couldn't find dir");
        return viqindex;
    }

    bool hasindex(const IQIndex& i) const { return findindex(i) != 0; }
    bool is_complex() const { return findindex(IQIndReIm) != 0; }
    bool has_virtual() const { return viqindex.index(1) != IQEmptyV.index(1); }
    QN virtualQN() const { return viqindex.qn(1); }
    const IQIndex& virtual_ind() const { return viqindex; }

    void addindex1(const IQIndex& I)
    {
        if(I.m() != 1) Error("IQTensor::operator*=(IQIndex): IQIndex must have m == 1.");    
        if(I.type() == Virtual) 
        {
            Arrow dir = I.dir();
            QN newq = I.qn(1);
            if(viqindex != IQEmptyV) //Add quantum numbers (with arrows)
            {
                newq = I.dir()*(viqindex.dir()*viqindex.qn(1)+I.dir()*I.qn(1));
                if(newq.Nf() < 0) //prefer to have Nf >= 0
                {
                    newq *= -1;
                    dir = dir*Switch;
                }
            }
            viqindex = IQIndex(I.name(),Index(I.index(1)),newq,dir);
            return;
        }
        solo(); p->uninit_rmap();
        foreach(ITensor& t, p->itensor) t.addindex1(I.index(1));
        p->iqindex_.push_back(I);
    }

    //----------------------------------------------------
    //IQTensor miscellaneous methods

    Real unique_Real() const
    {
        Real ur = 0.0;
        foreach(const IQIndex& I, p->iqindex_) ur += I.unique_Real();
        return ur;
    }
    int num_index(IndexType t) const 
    { 
        int count = 0;
        foreach(const IQIndex& I, p->iqindex_)
        if(I.type() == t) ++count;
        return count;
    }

    Real norm() const
    {
        Real res;
        foreach(const ITensor& t, p->itensor) res += sqr(t.norm());
        return sqrt(res);
    }

    int vec_size() const
    {
        int s = 0;
        for(const_iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
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
        solo();
        if(vec_size() != v.Length())
            Error("bad size");
        int off = 1;
        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
            {
            int d = jj->dat().Length();
            jj->AssignFromVec(v.SubVector(off,off+d-1));
            off += d;
            }
    }
    void GetSingComplex(Real& re, Real& im) const;

    
    void Randomize() { solo(); foreach(ITensor& t, p->itensor) t.Randomize(); }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    void printIQInds(string name = "") const
    { 
        cerr << "\nIQTensor ------------------\n";
        for(size_t j = 0; j < p->iqindex_.size(); ++j)
        { cerr << p->iqindex_[j] << "\n\n"; }
        if(has_virtual()) cerr << "Virtual = " << viqindex << "\n\n";
        cerr << "---------------------------\n\n";
    }

    void Assign(const IQTensor& other) const
    {
        map<ApproxReal,iten_it> semap;
        for(iten_it i = p->itensor.begin(); i != p->itensor.end(); ++i)
        { semap[ApproxReal(i->unique_Real())] = i; }
        for(const_iten_it i = other.p->itensor.begin(); i != other.p->itensor.end(); ++i)
        {
            ApproxReal se = ApproxReal(i->unique_Real());
            if(semap.count(se) == 0)
            {
                cout << "warning Assign semap.count is 0" << endl;
                cerr << "offending ITensor is " << *i << "\n";
                Error("bad Assign count se");
            }
            else { semap[se]->Assign(*i); }
        }
    }

    void SplitReIm(IQTensor& re, IQTensor& im) const;

    void conj()
    {
        solo();
        if(!is_complex())
        { foreach(IQIndex& I,p->iqindex_) I.conj(); }
        else
        {
            IQTensor r,i;
            SplitReIm(r,i);
            foreach(IQIndex& I,r.p->iqindex_) I.conj();
            foreach(IQIndex& I,i.p->iqindex_) I.conj();
            i *= -1.0;
            *this = r * IQComplex_1 + IQComplex_i * i;
        }
    }

private:
    void match_order() const
    {
        const int s = p->iqindex_.size();
        foreach(ITensor& t, p->itensor)
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
                if(p->iqindex_[k].hasindex(t.index(j)))
                { P.from_to(j,k+1); gotone = true; break; }
                if(!gotone) Error("match_order: !gotone");
            }
            t.Reshape(P);
        }
    }

public:

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

//Checks if the divergence of this IQTensor is zero
inline bool check_QNs(const ITensor& t) { return true; }
inline bool check_QNs(const IQTensor& T)
{
    foreach(const ITensor& it, T.itensors())
    {
        QN qtot;
        for(int j = 1; j <= it.r(); ++j) 
        { qtot += T.qn(it.index(j))*T.dir(it.index(j)); }
        qtot += T.virtual_ind().qn(1)*T.virtual_ind().dir();
        if(qtot != QN()) 
        {
            cerr << "check_QNs: IQTensor failed to have zero divergence.\n";
            cerr << "\nqtot = " << qtot << "\n\n";
            printdat = false; cerr << "Offending ITensor = " << it << "\n\n";
            cerr << "IQIndices = \n";
            foreach(const IQIndex& I, T.iqinds())
            { cerr << "\n" << I << "\n"; }
            return false;
        }
    }
    return true;
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

IQTensor IQTSing(makeSing);
IQTensor IQComplex_1(makeComplex_1), IQComplex_i(makeComplex_i);
const IQIndex& IQTensor::ReImIndex = IQIndReIm;
#endif
#endif
