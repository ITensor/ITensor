#ifndef __TENSOR_H
#define __TENSOR_H
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <assert.h>
#include <error.h> //utilities
#include <math.h>
#include <cstdlib>
#include <vector>
#include <list>
#include "boost/foreach.hpp"
#include "boost/format.hpp"
#include "boost/intrusive_ptr.hpp"
#include "boost/noncopyable.hpp"
//#include <tr1/array>
#include "boost/array.hpp"
#include "boost/uuid/uuid.hpp"
#include "boost/uuid/random_generator.hpp"
#include "boost/uuid/string_generator.hpp"
#include "boost/static_assert.hpp"
//#include "boost/phoenix/core.hpp"
//#include "boost/phoenix/operator.hpp"
//#include "boost/phoenix/core/reference.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::string;
using std::stringstream;
using std::ostringstream;
using std::setprecision;
using std::pair;
using std::make_pair;
using std::map;
using std::vector;
using std::list;
using boost::noncopyable;
using boost::intrusive_ptr;
using boost::format;
//using std::tr1::array;
using boost::array;
using namespace std::rel_ops;
using boost::uuids::uuid;
using boost::uuids::random_generator;
using boost::uuids::string_generator;
//using boost::phoenix::arg_names::arg1;
//using boost::phoenix::arg_names::arg2;
//using boost::phoenix::arg_names::arg3;
//using boost::phoenix::arg_names::arg4;
//using boost::phoenix::arg_names::arg5;
//using boost::phoenix::arg_names::arg6;
//using boost::phoenix::arg_names::arg7;
//using boost::phoenix::arg_names::arg8;
//using boost::phoenix::arg_names::arg9;
//using boost::phoenix::arg_names::arg10;
//using boost::phoenix::ref;

static const Real MAX_CUT = 1E-15;
static const int MAX_M = 5000;

Real ran1();

//#define COUNT_COPIES
//#define COLLECT_PRODSTATS

//----------------------------------
//For bounds checking - can remove once implementation is tested
//#define ITENSOR_USE_AT

#ifdef  ITENSOR_USE_AT
#define GET(container,j) (container.at(j))
#else
#define GET(container,j) (container[j])
#endif

//----------------------------------

#ifndef DNDEBUG
#define DO_IF_DEBUG(x) { x }
#else
#define DO_IF_DEBUG(x) { }
#endif

#ifdef COUNT_COPIES
#define IF_COUNT_COPIES(x) { x }
#else
#define IF_COUNT_COPIES(x) { }
#endif


#ifdef COLLECT_PRODSTATS
#include <cputime.h>
#endif

#define foreach BOOST_FOREACH

//---------------------------------------
#define ENABLE_INTRUSIVE_PTR(ClassName) \
friend inline void intrusive_ptr_add_ref(ClassName* p) { ++(p->numref); } \
friend inline void intrusive_ptr_release(ClassName* p) { if(--(p->numref) == 0){ delete p; } } \
int count() const { return numref; }
//---------------------------------------

const Real Pi = M_PI;
const Real Sqrt2 = sqrt(2);
const Real ISqrt2 = 1.0/sqrt(2);
const Real Sqrt3 = sqrt(3);
const Real ISqrt3 = 1.0/sqrt(3);

inline Real sqr(Real x) { return x*x; }

#define NMAX 8

extern bool printdat;
extern Vector lastd;
struct DMRGOpts; extern const DMRGOpts DefaultDMRGOpts;
extern bool catch_debug;
#ifdef COUNT_COPIES
extern int copycount;
#endif
#ifdef COLLECT_PRODSTATS
class Prodstats; extern Prodstats prodstats;
#endif

enum Direction { Fromright, Fromleft, Both, None };

enum Printdat { ShowData, HideData };

#define Print(X) { printdat = false; cerr << "\n" << #X << " =\n" << X << "\n"; }
#define PrintDat(X) { printdat = true; cerr << "\n" << #X << " =\n" << X << "\n"; printdat = false; }

//Enum defining directions for arrows
enum Arrow { In = -1, Out = 1 };

inline Arrow operator*(const Arrow& a, const Arrow& b)
{ return (int(a)*int(b) == In) ? In : Out; }

const Arrow Switch = In*Out;

inline ostream& operator<<(ostream& s, const Arrow& D)
{ if(D == In) s << "In"; else s << "Out"; return s; }


template<class T, class Op> void for_all(T& a, Op f) { for_each(a.begin(),a.end(),f); }

template<class T> vector<T>& operator*=(vector<T>& v1, const vector<T>& v2) 
{
    const unsigned int sz = v1.size();
    assert(v2.size() == sz);
    for(unsigned int n = 0; n < sz; ++n) v1[n] *= v2[n];
    return v1;
}
template<class T> vector<T> operator*(const vector<T>& v1, const vector<T>& v2) 
{ vector<T> res(v1); res *= v2; return res; }

template<class T> vector<T>& operator*=(vector<T>& v1, const vector<T*>& v2) 
{
    const unsigned int sz = v1.size();
    assert(v2.size() == sz);
    for(unsigned int n = 0; n < sz; ++n) v1[n] *= *(v2[n]);
    return v1;
}
template<class T> vector<T> operator*(const vector<T>& v1, const vector<T*>& v2) 
{ vector<T> res(v1); res *= v2; return res; }

template<class T> vector<T> operator*(const vector<const T*>& v1, const vector<const T*>& v2) 
{ 
    const unsigned int sz = v1.size();
    assert(v2.size() == sz);
    vector<T> res(sz); 
    for(unsigned int n = 0; n < sz; ++n) res[n] = *(v1[n]) * *(v2[n]);
    return res; 
}

template<class T>
ostream& operator<<(ostream& s, const vector<T>& v)
{ 
    if(v.size() == 0) s << "(Empty vector)\n";
    for(unsigned int n = 0; n < v.size(); ++n) { s << n << ": " << v[n] << "\n"; } 
    return s; 
}

template<class T> T& operator*=(T& t1, const T* pt2) 
{ t1 *= *(pt2); return t1; }
template<class T> T operator*(const T& t1, const T* pt2) 
{ T res(t1); res *= *(pt2); return res; }


class ApproxReal
{
public:
    Real r;
    ApproxReal(Real _r = 0.0) : r(_r) {}

    friend inline bool operator==(const ApproxReal &a,const ApproxReal &b)
    { return fabs(a.r-b.r) < 1.0e-12; }
    friend inline bool operator<(const ApproxReal &a,const ApproxReal &b)
    { return b.r-a.r > 1.0e-12; }
};


enum IndexType { Link, Site, ReIm, Virtual };
static const char * indextypename[] = { "Link","Site","ReIm","Virtual" };
enum PrimeType { primeLink, primeSite, primeBoth, primeNone };

inline ostream& operator<<(ostream& s, const IndexType& it)
{ 
    if(it == Link) s << "Link"; 
    else if(it == Site) s << "Site"; 
    else if(it == ReIm) s << "ReIm"; 
    else if(it == Virtual) s << "Virtual"; 
    return s; 
}

inline int IndexTypeToInt(IndexType it)
{
    if(it == Link) return 1;
    if(it == Site) return 2;
    if(it == ReIm) return 3;
    if(it == Virtual) return 4;
    Error("No integer value defined for IndexType.");
    return -1;
}
inline IndexType IntToIndexType(int i)
{
    if(i == 1) return Link;
    if(i == 2) return Site;
    if(i == 3) return ReIm;
    if(i == 4) return Virtual;
    cerr << format("No IndexType value defined for i=%d\n")%i,Error("");
    return Virtual;
}

inline string putprimes(string s, int plev = 0)
{ for(int i = 1; i <= plev; ++i) s += "\'"; return s;}

inline string nameindex(IndexType it, int plev = 0)
{ return putprimes(string(indextypename[(int)it]),plev); }

inline string nameint(string f,int ix)
{ stringstream ss; ss << f << ix; return ss.str(); }

enum Imaker {makeReIm,makeReImP,makeReImPP,makeEmptyV,makeNull};

#define UID_NUM_PRINT 2
inline ostream& operator<<(ostream& s, const uuid& id)
{ 
    s.width(2);
    for(uuid::size_type i = id.size()-UID_NUM_PRINT; i < id.size(); ++i) 
    {
        s << static_cast<unsigned int>(id.data[i]);
    }
    s.width(0);
    return s; 
}

struct UniqueID
{
    uuid id;

    UniqueID() : id(random_generator()()) { }

    UniqueID& operator++()
    {
        int i = id.size(); 
        while(--i >= 0)
        { 
            if(++id.data[i] == 0) continue; 
            break;
        }
        return *this;
    }

    operator uuid() const { return id; }

    friend inline ostream& operator<<(ostream& s, const UniqueID& uid) { s << uid.id; return s; }
};

inline int prime_number(int n)
{
    static const array<int,54> plist = { { 
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 
    79, 83, 89, 97, 101, 103, 107, 109, 113, 
    127, 131, 137, 139, 149, 151, 157, 163, 
    167, 173, 179, 181, 191, 193, 197, 199, 
    211, 223, 227, 229, 233, 239, 241, 251 
    } };
    return plist.at(n);
}


namespace Internal {

//Storage for Index's
class IndexDat : public noncopyable
{
    mutable unsigned int numref;
    const bool is_static_;
public:
    //static int indcount; //Global ID for IndexDats
    static UniqueID lastID;

    IndexType _type;
    //const int ind; 
    uuid ind;
    const int m_;
    Real ur;
    string sname;

    void set_unique_Real()
	{
        //ur = sin(ind * sqrt(1.0/7.0) + ((int)_type - (int)Site) * sqrt(1.0 / 13.0));
        Real arg = 0;
        int pn = 1;
        for(int i = int(ind.size()); i >= 0; --i)
        { arg += ind.data[i]*sqrt(1.0/(prime_number(++pn)*1.0)); }
        arg *= sqrt(1.0/(prime_number(++pn)*1.0));
        arg += ((int)_type - (int)Site) * sqrt(1.0/(prime_number(++pn)*1.0));
        ur = sin(arg);
        //cerr << format("Set a unique real of %.25f, hash_value(ind) = %d\n")%ur%hash_value(ind);
	}

    IndexDat(string name="", int mm = 1,IndexType it=Link) :
    numref(0), is_static_(false),
    _type(it), 
    //ind(++indcount), 
    ind(++lastID),
    m_(mm), 
    sname(name)
	{ 
        if(it == ReIm) Error("bad call to create IndexDat with type ReIm");
        assert((it==Virtual ? (mm==1) : true)); //If type is Virtual, m must be 1
        set_unique_Real();
	}

    //For use with read/write functionality of Index class
    IndexDat(string ss, int mm, IndexType it, uuid ind_) :
    numref(0), is_static_(false), _type(it), ind(ind_), m_(mm), sname(ss)
	{ 
        if(it == ReIm) Error("bad call to create IndexDat with type ReIm");
        assert((it==Virtual ? (mm==1) : true)); //If type is Virtual, m must be 1
        set_unique_Real();
	}

    //Don't actually use random uuid generator for these static IndexDats
    IndexDat(Imaker im) : 
    numref(1000000000), is_static_(true),
    _type(ReIm), 
    //ind((im==makeNull || im==makeEmptyV) ? 0 : 1), 
    m_( (im==makeNull || im==makeEmptyV) ? 1 : 2)
	{ 
        string_generator gen;
        if(im==makeNull || im==makeEmptyV) 
        { ind = gen("{00000000-0000-0000-0000-000000000000}"); }
        else                               
        { ind = gen("{10000000-0000-0000-0000-000000000000}"); }

        if(im == makeNull)
        {
            _type = Site;
            ur = 0.0;
            return;
        }
        else if(im == makeReIm) sname = "ReIm";
        else if(im == makeReImP) sname = "ReImP";
        else if(im == makeReImPP) sname = "ReImPP";
        else if(im == makeEmptyV) 
        {
            _type = Virtual;
            sname = "EmptyVirtual";
        }
        set_unique_Real(); 
	}

    friend inline void intrusive_ptr_add_ref(IndexDat* p) { ++(p->numref); }
    friend inline void intrusive_ptr_release(IndexDat* p) { if(!p->is_static_ && --(p->numref) == 0){ delete p; } }
    int count() const { return numref; }
};
} // namespace Internal

extern Internal::IndexDat IndexDatNull, IndReDat, IndReDatP, IndReDatPP, IndEmptyVDat;

struct IndexVal;

class Index
{
protected:
    intrusive_ptr<Internal::IndexDat> p;
public:
    int primelevel; 

    int m() const { return p->m_; }
    inline string showm() const { return (format("m=%d")%(p->m_)).str(); }
    //int Ind() const { return p->ind; }
    boost::uuids::uuid Ind() const { return p->ind; }
    IndexType type() const { return p->_type; }
    void settype(IndexType t) { p->_type = t; }
    string rawname() const { return p->sname; }
    Real unique_Real() const { assert(p!=0); return p->ur*(1+primelevel); }
    string name() const  { return putprimes(rawname(),primelevel); }
    void setname(string newname) { p->sname = newname; }
    bool is_null() const { return (p == &IndexDatNull); }
    int count() const { return p->count(); }

    Index() : p(&IndexDatNull), primelevel(0) { }

    Index(string name, int mm = 1, IndexType it=Link, int plev = 0) 
	: p(new Internal::IndexDat(name,mm,it)), primelevel(plev) { }

    Index(istream& s) { read(s); }


    Index(Imaker im)
	{
        if(im == makeNull)
            p = &IndexDatNull, primelevel = 0;
        else if(im == makeReIm)
            p = &IndReDat, primelevel = 0;
        else if(im == makeReImP)
            p = &IndReDatP,  primelevel = 1;
        else if(im == makeReImPP)
            p = &IndReDatPP,  primelevel = 2;
        else if(im == makeEmptyV)
            p = &IndEmptyVDat, primelevel = 0;
        else Error("Unrecognized Imaker type.");
	}

    Index(PrimeType pt,const Index& other, int primeinc = 1) 
	: p(other.p), primelevel(other.primelevel)
	{
        primelevel = other.primelevel;
        for(int i = 1; i <= primeinc; ++i) doprime(pt);
	}

    // rel_ops defines the other comparisons based on == and <
    bool operator==(const Index& other) const 
	{ return (unique_Real() == other.unique_Real()); }

    bool operator<(const Index& other) const 
	{ return (unique_Real() < other.unique_Real()); }

    IndexVal operator()(int i) const;

    bool noprime_equals(const Index& other) const
	{ return (p->_type == other.p->_type && p->ind == other.p->ind); }

    void mapprime(int plevold, int plevnew, PrimeType pr = primeBoth)
	{
        if(type() == ReIm) return;
        if(primelevel != plevold) return;
        else if( (pr == primeBoth && type() != Virtual)
        || (type() == Site && pr == primeSite) 
        || (type() == Link && pr == primeLink) )
        {
            primelevel = plevnew;
        }
	}
    void doprime(PrimeType pr, int inc = 1)
	{
        if(type() == ReIm) return;
        if( (pr == primeBoth && type() != Virtual)
        || (type() == Site && pr == primeSite) 
        || (type() == Link && pr == primeLink) )
        {
            primelevel += inc;
        }
	}
    Index primed(int inc = 1) const { return Index(primeBoth,*this,inc); }

    Index deprimed() const { Index cp(*this); cp.primelevel = 0; return cp; }

    void noprime(PrimeType p = primeBoth) { doprime(p,-primelevel); }

    friend inline ostream & operator << (ostream & s, const Index & t)
    {
        if(t.name() != "" && t.name() != " ") s << t.name() << "/";
        return s << nameindex(t.type(),t.primelevel) << "-" << t.Ind() << ":" << t.m();
    }

    void write(ostream& s) const 
    { 
        if(is_null()) Error("Index::write: Index is null");
        s.write((char*) &primelevel,sizeof(primelevel));
        const int t = IndexTypeToInt(p->_type);
        s.write((char*) &t,sizeof(t));
        //s.write((char*) &(p->ind),sizeof(p->ind));
        for(int i = 0; i < int(p->ind.size()); ++i) 
        { const char c = p->ind.data[i] - '0'; s.write(&c,sizeof(c)); }
        s.write((char*) &(p->m_),sizeof(p->m_));
        const int nlength = p->sname.length();
        s.write((char*) &nlength,sizeof(nlength));
        s.write(p->sname.data(),nlength+1);
    }

    void read(istream& s)
    {
        s.read((char*) &primelevel,sizeof(primelevel));
        int t; s.read((char*) &t,sizeof(t));
        //int ind; s.read((char*) &ind,sizeof(ind));
        boost::uuids::uuid ind;
        for(int i = 0; i < int(ind.size()); ++i) 
        { char c; s.read(&c,sizeof(c)); ind.data[i] = '0'+c; }
        int mm; s.read((char*) &mm,sizeof(mm));
        int nlength; s.read((char*) &nlength,sizeof(nlength));
        char* newname = new char[nlength+1]; s.read(newname,nlength+1);
        string ss(newname); delete newname;
        p = new Internal::IndexDat(ss,mm,IntToIndexType(t),ind);
    }

    void print(string name = "") const
    { cerr << "\n" << name << " =\n" << *this << "\n"; }

    void conj() { } //for forward compatibility with arrows

}; //class Index
extern Index IndNull, IndReIm, IndReImP, IndReImPP, IndEmptyV;

template <class T> 
T conj(T res) { res.conj(); return res; }

class ITensor;

struct IndexVal
{
    Index ind; 
    int i;
    IndexVal() : ind(IndNull),i(0) { }
    IndexVal(const Index& index, int i_) : ind(index),i(i_) { assert(i <= ind.m()); }
    inline friend ostream& operator<<(ostream& s, const IndexVal& iv)
    { return s << "IndexVal: i = " << iv.i << ", ind = " << iv.ind << "\n"; }
    ITensor operator*(const IndexVal& oth) const;
    ITensor operator*(Real fac) const;
    friend inline ITensor operator*(Real fac, const IndexVal& iv);
    IndexVal primed() const { return IndexVal(ind.primed(),i); }
};
extern IndexVal IVNull;

class Permutation // Tell where each index will go, p(2,1,3) says 1 -> 2, 2 -> 1, 3 -> 3
{
    typedef array<int,NMAX+1> int9;
    void set8(int9 *n, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8)
    {
        (*n)[1] = i1; (*n)[2] = i2; (*n)[3] = i3; (*n)[4] = i4;
        (*n)[5] = i5; (*n)[6] = i6; (*n)[7] = i7; (*n)[8] = i8;
    }
public:
    array<int,NMAX+1> ind;
    Permutation(int i1 = 1, int i2 = 2, int i3 = 3, int i4 = 4, int i5 = 5, int i6 = 6,
	    int i7 = 7,  int i8 = 8)
	{ set8(&ind,i1,i2,i3,i4,i5,i6,i7,i8); }
    void check(int d)
	{
        for(int i = 1; i <= d; i++)
        if(ind[i] > d || ind[i] < 1) Error("bad Permutation level 1");

        for(int i = 1; i <= d; i++)
        for(int j = 1; j <= d; j++)
        if(ind[i] == ind[j] && i != j) Error("bad Permutation level 2");
	}

    Permutation inverse() const
    {
        Permutation inv;
        for(int n = 1; n <= NMAX; ++n) inv.ind[ind[n]] = n; 
        return inv;
    }

    friend inline ostream& operator<<(ostream& s, const Permutation& p)
    {
        for(int i = 1; i <= NMAX; i++) s << format("(%d,%d) ") % i % p.ind[i];
        return s;
    }
};


enum ITmaker {makeComplex_1,makeComplex_i,makeConjTensor};

namespace Internal {

//#define DO_ALT

#ifdef DO_ALT
struct PDat
{
    Permutation I; 
    Vector v;
    PDat(const Permutation& P_, const Vector& v_) : I(P_.inverse()), v(v_) { }
    PDat(const Permutation& P_) : I(P_.inverse()) { }
};
#endif

//Storage for ITensors
class ITDat : noncopyable
{
private:
    mutable unsigned int numref;
public:
    Vector v;
#ifdef DO_ALT
    vector<PDat> alt;
#endif

    ITDat(int size) : numref(0), v(size)
	{ if(size > 0) v = 0; }

    ITDat(istream& s) : numref(0)
	{ 
        int size = 0;
        s.read((char*) &size,sizeof(size));
        v.ReDimension(size);
        s.read((char*) v.Store(), sizeof(Real)*size);
    }

    void write(ostream& s) const 
    { 
        const int size = v.Length();
        s.write((char*) &size, sizeof(size));
        s.write((char*) v.Store(), sizeof(Real)*size); 
    }

    void print() const { cout << "ITDat: v = " << v; }

    friend class ITensor;
    ENABLE_INTRUSIVE_PTR(ITDat)
private:
    ~ITDat() { } //must be dynamically allocated
};

} //namespace Internal

class ITensor //Index Tensor
{
private:
    mutable intrusive_ptr<Internal::ITDat> p; //mutable: const methods may want to reshape data
    int rn;
    mutable array<Index,NMAX+1> _indexn; //Indices having m!=1, maximum of 8 (_indexn[0] not used), mutable to allow reordering
    mutable vector<Index>       _index1; //Indices having m==1
    mutable Real _logfac; //mutable since e.g. normlogto is logically const
    Real ur;
    mutable bool _neg; //true if overall sign is -1, mutable since e.g. solo_dosign logically const

    void allocate(int dim) { p = new Internal::ITDat(dim); }

#ifdef DO_ALT
    void newAltDat(const Permutation& P) const
    { p->alt.push_back(PDat(P)); }
    PDat& lastAlt() const { return p->alt.back(); } 
#endif

    //Disattach self from current ITDat and create own copy instead.
    //Necessary because ITensors logically represent distinct
    //objects even though they may share data in reality.
    void solo_dosign() const
	{
        assert(p != 0);
        if(p->count() != 1) 
        {
            intrusive_ptr<Internal::ITDat> new_p = new Internal::ITDat(0);
            new_p->v = p->v;
            p.swap(new_p);
            IF_COUNT_COPIES(++copycount;)
        }
        if(_neg) 
        { 
            p->v *= -1; 
#ifdef DO_ALT
            foreach(PDat& pd, p->alt) pd.v *= -1;
#endif
            _neg = false; 
        }
	}

    void set_dat(const Vector& newv)
	{
        assert(p != 0);
        if(p->count() != 1) { intrusive_ptr<Internal::ITDat> new_p = new Internal::ITDat(0); p.swap(new_p); }
        p->v = newv;
#ifdef DO_ALT
        p->alt.clear();
#endif
	}
    
    void set_unique_Real()
	{
        ur = 0.0;
        for(int j = 1; j <= rn; ++j)
        { ur += GET(_indexn,j).unique_Real(); }
        foreach(const Index& I, _index1)
        { ur += I.unique_Real(); }
	}

    void _construct2(const Index& i1, const Index& i2)
    {
        if(i1.m()==1) _index1.push_back(i1); else { GET(_indexn,++rn) = i1; }
        if(i2.m()==1) _index1.push_back(i2); else { GET(_indexn,++rn) = i2; }
        allocate(i1.m()*i2.m()); 
        set_unique_Real();
    }

    //Prefer to map via a Combiner
    //'mapindex' is useful and efficient if used in a safe way, however.
    void mapindex(const Index& i1, const Index& i2)
	{
        assert(i1.m() == i2.m());
        solo_dosign();
        for(int j = 1; j <= rn; ++j) 
        if(GET(_indexn,j) == i1) 
        {
            GET(_indexn,j) = i2;
            set_unique_Real();
            return;
        }
        vector<Index>::iterator i1pos = find(_index1.begin(),_index1.end(),i1);
        if(i1pos == _index1.end()) 
        {
            cerr << "\nFor ITensor = " << *this << "\n";
            cerr << "Missing index i1 = " << i1 << "\n";
            Error("ITensor::mapindex(i1,i2): ITensor does not have index i1");
        }
        *i1pos = i2;
        set_unique_Real();
	}

    void getperm(const ITensor& other, Permutation& P) const;

    friend void toMatrixProd(const ITensor& L, const ITensor& R, 
                             array<bool,NMAX+1>& contractedL, array<bool,NMAX+1>& contractedR, 
                             MatrixRefNoLink& lref, MatrixRefNoLink& rref);
public:
    typedef Index IndexT;
    typedef IndexVal IndexValT;

    //Accessor Methods ----------------------------------------------

    Real unique_Real() const { return ur; }	// depends on indices only, unordered
    const Index& index(int j) const { return (j > rn ? GET(_index1,j-(rn+1)) : GET(_indexn,j)); }
    const Index& indexn(int j) const { return GET(_indexn,j); }
    int r() const { return rn + _index1.size(); }
    int r_n() const { return rn; }
    int r_1() const { return _index1.size(); }
    int m(int j) const { return (j > rn ? 1 : _indexn[j].m()); }

    bool is_null() const { return (p == 0); }
    bool is_not_null() const { return (p != 0); }
    bool is_complex() const { return findindexn(IndReIm) > 0; }
    bool is_not_complex() const { return (findindexn(IndReIm) == 0); }
    Vector& ncdat() { assert(p != 0); solo_dosign(); return p->v; } //Can we get dat & ncdat to do the right thing automatically?
    const VectorRef dat() const { assert(p != 0); return p->v*(_neg ? -1 : 1); }
    int Length() const { return dat().Length(); }
    Real logfac() const { return _logfac; }
    bool neg() const { return _neg; }
    void setlogfac(Real newlogfac) { _logfac = newlogfac; }

    //These methods can be used for const iteration over Indices in a foreach loop
    //e.g. foreach(const Index& I, t.indexn() ) { ... }
    typedef array<Index,NMAX+1>::const_iterator indexn_it;
    const pair<indexn_it,indexn_it> indexn() const { return make_pair(_indexn.begin()+1,_indexn.begin()+rn+1); }
    const vector<Index>&            index1() const { return _index1; }

    //Constructors --------------------------------------------------

    ITensor() : p(0), rn(0), _logfac(0), ur(0), _neg(false)  { }

    ITensor(istream& s) { read(s); }

    ITensor(Real val) : rn(0), _logfac(0), _neg(false)
	{ 
        allocate(1);
        p->v = val;
        set_unique_Real();
    }

    ITensor(const Index& i1) : rn(0), _logfac(0), _neg(false)
	{ 
        if(i1.m()==1) _index1.push_back(i1); else { _indexn[1] = i1; ++rn; }
        allocate(i1.m());
        set_unique_Real();
    }

    ITensor(const Index& i1, Real val) : rn(0), _logfac(0), _neg(false)
	{ 
        if(i1.m()==1) _index1.push_back(i1); else { _indexn[1] = i1; ++rn; }
        allocate(i1.m()); p->v = val; 
        set_unique_Real();
    }

    ITensor(const Index& i1, const Vector& V) : rn(0), _logfac(0), _neg(false)
	{ 
        if(i1.m() != V.Length()) Error("Mismatch of Index and Vector sizes.");
        if(i1.m()==1) _index1.push_back(i1); else { GET(_indexn,1) = i1; ++rn; }
        allocate(0); p->v = V; 
        set_unique_Real();
    }

    ITensor(Index i1,Index i2) : rn(0), _logfac(0), _neg(false)
	{ _construct2(i1,i2); }

    //Create an ITensor as a matrix with 'a' on the diagonal
    ITensor(Index i1,Index i2,Real a) : rn(0), _logfac(0), _neg(false)
    {
        _construct2(i1,i2);
        int nn = min(i1.m(),i2.m());
        if(rn == 2) for(int i = 1; i <= nn; ++i) p->v((i-1)*i1.m()+i) = a;
        else p->v(1) = a;
    }

    ITensor(Index i1,Index i2,const MatrixRef& M) : rn(0), _logfac(0), _neg(false)
    {
        _construct2(i1,i2);
        if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
        { Error("ITensor(Index,Index,Matrix): Mismatch of Index sizes and matrix."); }
        MatrixRef dref; p->v.TreatAsMatrix(dref,i2.m(),i1.m());
        dref = M.t();
    }

    ITensor(Index i1, Index i2, Index i3,
            Index i4 = IndNull, Index i5 = IndNull, Index i6 = IndNull,
            Index i7 = IndNull, Index i8 = IndNull)
            : rn(0), _logfac(0), _neg(false)
    {
        array<Index,NMAX+1> ii = {{ IndNull, i1, i2, i3, i4, i5, i6, i7, i8 }};
        int dim = 1;
        for(int n = 1; n <= NMAX; ++n)
        { 
            if(ii[n] == IndNull) break;
            if(ii[n].m()==1) _index1.push_back(ii[n]); else { dim *= ii[n].m(); _indexn[++rn] = ii[n]; } 
        }
        allocate(dim);
        set_unique_Real();
    }

    ITensor(const IndexVal& iv, Real fac = 1) : rn(0), _logfac(0), _neg(false)
    { 
        if(iv.ind.m()==1) _index1.push_back(iv.ind); else { _indexn[++rn] = iv.ind; } 
        allocate(iv.ind.m());  
        p->v(iv.i) = fac; 
        set_unique_Real(); 
    }

    ITensor(IndexVal iv1, IndexVal iv2, IndexVal iv3 = IVNull,
            IndexVal iv4 = IVNull, IndexVal iv5 = IVNull, IndexVal iv6 = IVNull,
            IndexVal iv7 = IVNull, IndexVal iv8 = IVNull)
            : rn(0), _logfac(0), _neg(false)
	{
        array<IndexVal,NMAX+1> iv = {{ IVNull, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};
        int dim = 1;
        for(int n = 1; n <= NMAX; ++n)
        { 
            if(iv[n].ind == IndNull) break;
            if(iv[n].ind.m()==1) _index1.push_back(iv[n].ind); else { dim *= iv[n].ind.m(); _indexn[++rn] = iv[n].ind; } 
        }
        allocate(dim);
        //p->v((((((((iv8.i-1)*m(7)+iv7.i-1)*m(6)+iv6.i-1)*m(5)+iv5.i-1)*m(4)+iv4.i-1)*m(3)+iv3.i-1)*m(2)+iv2.i-1)*m(1)+iv1.i)
        //= 1;
        operator()(iv1.i,iv2.i,iv3.i,iv4.i,iv5.i,iv6.i,iv7.i,iv8.i) = 1;
        set_unique_Real();
    }

    template <typename IndexContainer>
    ITensor(const IndexContainer& I, bool do_allocate) : rn(0), _logfac(0), _neg(false)
    {
        int alloc_size = 1;
        foreach(const Index& i, I)
        {
            if(i == IndNull) Error("ITensor: null Index in constructor.");
            if(i.m()==1) _index1.push_back(i); 
            else 
            { 
                if(rn == NMAX) Error("ITensor(const vector<Index>& I): too many indices with m > 1");
                GET(_indexn,++rn) = i; 
                alloc_size *= i.m(); 
            }
        }
        if(do_allocate) allocate(alloc_size);
        else allocate(0);

        set_unique_Real();
    }

    ITensor(PrimeType pt, ITensor other) : rn(other.rn),
	_indexn(other._indexn), _index1(other._index1), _logfac(other._logfac), _neg(other._neg)
	{
        p.swap(other.p); //copy is made on passing the arg
        doprime(pt); 
        set_unique_Real();
	}
	
    ITensor(ITmaker itm) : rn(1), _logfac(0), _neg(false)
	{
        GET(_indexn,1) = IndReIm; allocate(2);
        if(itm == makeComplex_1)  { p->v(1) = 1; }
        if(itm == makeComplex_i)  { p->v(2) = 1; }
        if(itm == makeConjTensor) { p->v(1) = 1; p->v(2) = -1; }
        set_unique_Real();
	}

    //ITensor(const ITensor& other) - needed? See 'rule of two' article

    //Operators -------------------------------------------------------

    ITensor& operator*=(const ITensor& other);
    ITensor operator*(ITensor other) const { other *= *this; return other; }

    ITensor& operator*=(const IndexVal& iv) { ITensor oth(iv); return operator*=(oth); } 
    ITensor operator*(const IndexVal& iv) const { ITensor res(*this); res *= iv; return res; }
    friend inline ITensor operator*(const IndexVal& iv, ITensor t) { return (t *= iv); }

    ITensor& operator*=(Real fac) { _neg ^= (fac < 0); _logfac += log(fabs(fac)+1E-100); return *this; }
    ITensor operator*(Real fac) const { ITensor res(*this); res *= fac; return res; }
    friend inline ITensor operator*(Real fac, ITensor t) { return (t *= fac); }

    //operator '/' is actually non-contracting product
    ITensor& operator/=(const ITensor& other);
    ITensor operator/(const ITensor& other) const { ITensor res(*this); res /= other; return res; }

    ITensor& operator+=(const ITensor& o);
    ITensor operator+(const ITensor& o) const { ITensor res(*this); res += o; return res; }

    ITensor& operator-=(const ITensor& o)
    {
        if(this == &o) { _logfac -= 200; return *this; }
        ncdat() *= -1; operator+=(o); p->v *= -1; return *this; 
    }
    ITensor operator-(const ITensor& o) const { ITensor res(*this); res -= o; return res; }

    //ITensor& operator=(const ITensor& other) - needed? See 'rule of two' article

    //Index Methods ---------------------------------------------------

    Index findtype(IndexType t) const
	{
        for(int j = 1; j <= rn; ++j)
        if(GET(_indexn,j).type() == t) return GET(_indexn,j);
        foreach(const Index& I,_index1)
        if(I.type() == t) return I;
        Error("ITensor::findtype failed."); return Index();
	}

    bool findtype(IndexType t,Index& I) const
	{
        for(int j = 1; j <= rn; ++j)
        if(GET(_indexn,j).type() == t)
        {
            I = GET(_indexn,j);
            return true;
        }
        foreach(const Index& J,_index1)
        if(J.type() == t)
        {
            I = J;
            return true;
        }
        return false;
	}

    int findindex(const Index& I) const
    {
        if(I.m() == 1) return findindex1(I);
        else           return findindexn(I);
        return 0;
    }

    int findindexn(const Index& I) const
	{
        if(I.m() == 1) return 0;
        for(int j = 1; j <= rn; ++j)
        if(GET(_indexn,j) == I) return j;
        return 0;
	}

    int findindex1(const Index& I) const
	{
        if(I.m() != 1) return 0;
        for(unsigned int j = 0; j < _index1.size(); ++j)
        if(_index1[j] == I) return j+rn+1;
        return 0;
	}

    //Checks that if this has an Index I, then it
    //also has I' (I' can have any primelevel)
    bool has_symmetric_nindices() const
    {
        for(int i = 1; i <= rn; ++i)
        for(int j = 1; j <= rn; ++j)
        {
            if(j == i) continue;
            if(GET(_indexn,i).noprime_equals(GET(_indexn,j))) break;
            if(j == rn) return false;
        }
        return true;
    }
    
    bool hasindex(const Index& I) const
	{
        if(I.m() == 1) return hasindex1(I);
        else           return hasindexn(I);
        return false;
	}

    bool hasindexn(const Index& I) const
	{
        for(int j = 1; j <= rn; ++j)
        if(GET(_indexn,j) == I) return true;
        return false;
	}

    bool hasindex1(const Index& I) const
	{
        foreach(const Index& J,_index1)
        if(J == I) return true;
        return false;
	}

    bool notin(const Index& I) const { return (findindexn(I)==0 && !hasindex1(I)); }

    template <class Iterable>
    void addindex1(const Iterable& indices) 
    { 
        if(_index1.empty()) { _index1.insert(_index1.begin(),indices.begin(),indices.end()); }
        else
        {
            _index1.insert(_index1.begin(),indices.begin(),indices.end());
            sort(_index1.begin(),_index1.end());
            vector<Index>::iterator it = unique(_index1.begin(),_index1.end());
            _index1.resize(it-_index1.begin());
        }
        set_unique_Real();
    }

    void addindex1(const Index& I) 
    { 
        if(I.m() != 1) { cerr << "I = " << I << "\n"; Error("ITensor::addindex1: m != 1"); }
        if(hasindex1(I)) return;
        _index1.push_back(I); 
        set_unique_Real();
    }

    void removeindex1(const Index& I) 
    { 
        vector<Index>::iterator it = find(_index1.begin(),_index1.end(),I);
        if(it == _index1.end()) Error("Couldn't find m == 1 Index to remove.");
        _index1.erase(it);
        set_unique_Real();
    }

    //Removes the jth index as found by findindex
    void removeindex1(int j) 
    { 
        vector<Index>::iterator it = _index1.begin()+(j-(rn+1));
        _index1.erase(it);
        set_unique_Real();
    }


    //Primelevel Methods ------------------------------------

    void noprime(PrimeType p = primeBoth)
	{
        for(int j = 1; j <= rn; ++j) _indexn[j].noprime(p);
        foreach(Index& I, _index1) I.noprime(p);
        set_unique_Real();
	}

    void doprime(PrimeType pt, int inc = 1)
	{
        for(int j = 1; j <= rn; ++j) _indexn[j].doprime(pt,inc);
        foreach(Index& I, _index1) I.doprime(pt,inc);
        set_unique_Real();
	}

    void primeall() { doprime(primeBoth,1); }
    void primesite() { doprime(primeSite,1); }
    void primelink() { doprime(primeLink,1); }

    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
	{
        for(int j = 1; j <= rn; ++j) GET(_indexn,j).mapprime(plevold,plevnew,pt);
        foreach(Index& I, _index1) I.mapprime(plevold,plevnew,pt);
        set_unique_Real();
	}

    void mapprimeind(const Index& I, int plevold, int plevnew, PrimeType pt = primeBoth)
	{
        for(int j = 1; j <= rn; ++j) 
        if(I == _indexn[j])
        {
            _indexn[j].mapprime(plevold,plevnew,pt);
            set_unique_Real();
            return;
        }
        foreach(Index& J, _index1) 
        if(I == J)
        {
            J.mapprime(plevold,plevnew,pt);
            set_unique_Real();
            return;
        }
        Print(*this);
        Print(I);
        Error("ITensor::mapprimeind: index not found.");
	}

    void primeind(const Index& I, int inc = 1) { mapindex(I,I.primed(inc)); }

    void primeind(const Index& I, const Index& J)
	{ 
        mapindex(I,I.primed());
        mapindex(J,J.primed());
	}

    void noprimeind(const Index& I) { mapindex(I,I.deprimed()); }

    void PrimesToBack();

    friend inline ITensor primed(ITensor A)
    { A.doprime(primeBoth,1); return A; }

    friend inline ITensor primesite(ITensor A)
    { A.doprime(primeSite,1); return A; }

    friend inline ITensor primelink(const ITensor& A)
    { ITensor res(A); res.doprime(primeLink,1); return res; }

    friend inline ITensor primeind(ITensor A, const Index& I)
    { A.mapindex(I,I.primed()); return A; }

    friend inline ITensor primeind(ITensor A, const Index& I1, const Index& I2)
    { A.mapindex(I1,I1.primed()); A.mapindex(I2,I2.primed()); return A; }

    friend inline ITensor deprimed(ITensor A)
    { A.noprime(); return A; }

    //Element Access Methods ----------------------------------------

    //Doesn't put in logfac or sign (i.e. _neg)
    Real operator()(int i1 = 1,int i2 = 1,int i3 = 1,int i4 = 1,int i5 = 1,
	    int i6 = 1,int i7 = 1, int i8 = 1) const
	{ assert(p != 0); return p->v((((((((i8-1)*m(7)+i7-1)*m(6)+i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1); }

    Real& operator()(int i1 = 1,int i2 = 1,int i3 = 1,int i4 = 1,int i5 = 1,
	    int i6 = 1,int i7 = 1, int i8 = 1) 
	{ assert(p != 0); solo_dosign(); return p->v((((((((i8-1)*m(7)+i7-1)*m(6)+i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1); }

    /*
    Real& val7(int i1,int i2,int i3,int i4,int i5,int i6,int i7)
	{ assert(rn==7); assert(p != 0); solo_dosign(); return p->v(((((((i7-1)*m(6)+i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1); }

    Real& val6(int i1,int i2,int i3,int i4,int i5,int i6)
	{ assert(rn==6); assert(p != 0); solo_dosign(); return p->v((((((i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1); }

    Real& val5(int i1,int i2,int i3,int i4,int i5)
	{ assert(rn==5); assert(p != 0); solo_dosign(); return p->v(((((i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1); }

    Real& val4(int i1,int i2,int i3,int i4)
	{ assert(rn==4); assert(p != 0); solo_dosign(); return p->v((((i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1); }

    Real& val3(int i1,int i2,int i3)
	{ assert(rn==3); assert(p != 0); solo_dosign(); return p->v(((i3-1)*m(2)+i2-1)*m(1)+i1); }
    */

    Real& val0()
	{ assert(rn==0); assert(p != 0); solo_dosign(); return p->v(1); }

    Real& val1(int i1)
	{ assert(rn==1); assert(p != 0); solo_dosign(); return p->v(i1); }

    Real val1(int i1) const
	{ assert(rn==1); assert(p != 0); return p->v(i1); }

    Real& val2(int i1,int i2)
	{ assert(rn==2); assert(p != 0); solo_dosign(); return p->v((i2-1)*m(1)+i1); }

    Real& operator()(const IndexVal& iv1, const IndexVal& iv2 = IVNull, const IndexVal& iv3 = IVNull,
                     const IndexVal& iv4 = IVNull, const IndexVal& iv5 = IVNull, const IndexVal& iv6 = IVNull,
                     const IndexVal& iv7 = IVNull, const IndexVal& iv8 = IVNull)
	{
        array<IndexVal,NMAX+1> iv = {{ IVNull, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};
        vector<int> ja(NMAX+1,1);
        int numgot = 0;
        for(int k = 1; k <= rn; ++k) //loop over indices of this ITensor
        {
            for(int j = 1; j <= NMAX; ++j)  // loop over the given indices
            if(_indexn[k] == iv[j].ind) 
            { ++numgot; ja[k] = iv[j].i; break; }
        }
        if(numgot != rn) 
        {
            cerr << format("numgot = %d, rn = %d\n")%numgot%rn;
            Error("ITensor::operator(): Not enough indices");
        }
	    assert(p != 0); solo_dosign(); 
        normlogto(0);
        return p->v((((((((ja[8]-1)*m(7)+ja[7]-1)*m(6)+ja[6]-1)*m(5)+ja[5]-1)*m(4)+ja[4]-1)*m(3)+ja[3]-1)*m(2)+ja[2]-1)*m(1)+ja[1]);
	}

    //Methods for Mapping to Other Objects --------------------------------------

    void Assign(const ITensor& other); // Assume *this and other have same indices but different order.
    		// Copy other into *this, without changing the order of indices in either
    		// operator= would put the order of other into *this

    void toMatrix11(const Index& i1, const Index& i2, Matrix& res, Real& logfac) const; //doesn't put in logfac
    void toMatrix11(const Index& i1, const Index& i2, Matrix& res) const; //puts in logfac
    void fromMatrix11(const Index& i1, const Index& i2, const Matrix& res);

    // group i1,i2; i3,i4
    void toMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,Matrix& res, Real& logfac) const; //doesn't put in logfac
    void toMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,Matrix& res) const; //puts in logfac
    void fromMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,const Matrix& res);

    // group i1,i2; i3
    void toMatrix21(const Index& i1, const Index& i2, const Index& i3, Matrix& res, Real& logfac) const; //doesn't put in logfac
    void toMatrix21(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const; //puts in logfac
    void fromMatrix21(const Index& i1, const Index& i2, const Index& i3, const Matrix& res);

    // group i1; i2,i3
    void toMatrix12(const Index& i1, const Index& i2, const Index& i3, Matrix& res, Real& logfac) const; //doesn't put in logfac
    void toMatrix12(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const; //puts in logfac
    void fromMatrix12(const Index& i1, const Index& i2, const Index& i3, const Matrix& res);

    int vec_size() const { return Length(); }
    void AssignToVec(VectorRef v) const
	{
        if(Length() != v.Length()) Error("ITensor::AssignToVec bad size");
        v = dat();
        v *= (_neg ? -1 : 1)*exp(_logfac);
	}
    void AssignFromVec(const VectorRef& v)
	{
        if(Length() != v.Length()) Error("ITensor::AssignToVec bad size");
        _logfac = 0; _neg = false;
        set_dat(v);
	}
    void ReshapeDat(const Permutation& p, Vector& rdat) const;
    void Reshape(const Permutation& p, ITensor& res) const;

    void read(istream& s)
    { 
        bool null_;
        s.read((char*) &null_,sizeof(null_));
        if(null_) { *this = ITensor(); return; }
        s.read((char*) &rn,sizeof(rn));
        s.read((char*) &_logfac,sizeof(_logfac));
        s.read((char*) &_neg,sizeof(_neg));
        size_t i1size = 0;
        s.read((char*) &i1size,sizeof(i1size));
        p = new Internal::ITDat(s);
        for(int j = 1; j <= rn; ++j) _indexn[j] = Index(s);
        _index1.reserve(i1size); for(size_t i = 1; i <= i1size; ++i) _index1.push_back(Index(s));
        set_unique_Real();
    }

    void write(ostream& s) const 
    { 
        bool null_ = is_null();
        s.write((char*) &null_,sizeof(null_));
        if(null_) return;
        s.write((char*) &rn,sizeof(rn));
        s.write((char*) &_logfac,sizeof(_logfac));
        s.write((char*) &_neg,sizeof(_neg));
        const size_t i1size = _index1.size();
        s.write((char*) &i1size,sizeof(i1size));
        p->write(s);
        for(int j = 1; j <= rn; ++j) _indexn[j].write(s);
        foreach(const Index& I, _index1) I.write(s);
    }

    //Other Methods -------------------------------------------------

    void Randomize() { ncdat().Randomize(); }

    void SplitReIm(ITensor& re, ITensor& im) const;
    void conj();
    //friend inline ITensor conj(ITensor A) { A.conj(); return A; }

    inline bool is_zero() const { return (norm() < 1E-20); } 

    Real norm() const { return Norm(dat()) * exp(_logfac); }

    void normalize() {  operator*=(1.0/norm()); }

    Real lognorm() const { return log(Norm(dat()) + 1.0e-100) + _logfac; }

    void donormlog()
	{
        Real f = Norm(dat());
        if(f != 0) { ncdat() *= 1.0/f; _logfac += log(f); }
	}

    void normlogto(Real newlogfac) const
	{
        Real dellogfac = newlogfac - _logfac;
        assert(p != 0); solo_dosign();
        if(dellogfac > 100.) p->v = 0;
        else                 p->v *= exp(-dellogfac);
        _logfac = newlogfac;
	}

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    friend ostream & operator << (ostream & s, const ITensor & t);

    bool checkDim() const
    {
        int dim = 1;
        for(int j = 1; j <= rn; ++j) dim *= _indexn[j].m();

        foreach(const Index& I, _index1)
        if(I.m() != 1) { cerr << I << "\n"; cerr << "WARNING: m != 1 index in _index1\n"; return false; }

        if(dim != Length()) 
        {
            print();
            cerr << "WARNING: Mismatched dat Length and Index dimension.\n";
            return false;
        }
        return true;
    }

}; //ITensor

extern ITensor Complex_1, Complex_i, ConjTensor;

inline Real Dot(const ITensor& x, const ITensor& y, bool doconj = true)
{
    if(x.is_complex())
	{
        ITensor res = (doconj ? conj(x) : x); res *= y;
        if(res.r() != 1) Error("Bad Dot 234234");
        return res(1)*exp(res.logfac());
	}
    else if(y.is_complex())
	{
        ITensor res = x; res *= y;
        if(res.r() != 1) Error("Bad Dot 37298789");
        return res(1)*exp(res.logfac());
	}

    ITensor res = x; res *= y;
    if(res.r() != 0) 
	{ x.print("x"); y.print("y"); Error("bad Dot"); }
    return res(1)*exp(res.logfac());
}

inline void Dot(const ITensor& x, const ITensor& y, Real& re, Real& im, bool doconj = true)
{
    if(x.is_complex())
	{
        ITensor res = (doconj ? conj(x) : x); res *= y;
        if(res.r() != 1) error("Bad Dot 334234");
        re = res(IndReIm(1)) * exp(res.logfac());
        im = res(IndReIm(2)) * exp(res.logfac());
        return;
	}
    else if(y.is_complex())
	{
        ITensor res = x; res *= y;
        if(res.r() != 1) error("Bad Dot 47298789");
        re = res(IndReIm(1)) * exp(res.logfac());
        im = res(IndReIm(2)) * exp(res.logfac());
        return;
	}
    if(x.r() != y.r()) 
	{
        cerr << "x = " << x << "\n";
        cerr << "y = " << y << "\n";
        Error("bad Dot 122414");
	}
    ITensor res = x; res *= y;
    if(res.r() != 0) 
	{
        cerr << "x = " << x << "\n";
        cerr << "y = " << y << "\n";
        Error("bad Dot 20234");
	}
    re = res(1)*exp(res.logfac());
    im = 0;
}

Index index_in_common(const ITensor& A, const ITensor& B, IndexType t);

class Counter
{
private:
    void init(int a)
	{
        i[0]=i[1]=i[2]=i[3]=i[4]=i[5]=i[6]=i[7]=i[8]=a;
        ind = 1;
	}
public:
    array<int,NMAX+1> n;
    array<int,NMAX+1> i;
    int r, ind;

    Counter() : r(0)
	{
        n[0] = 0;
        n[1]=n[2]=n[3]=n[4]=n[5]=n[6]=n[7]=n[8]=1;
        init(0);
	}

    Counter(const ITensor& t) : r(t.r())
	{
        n[0] = 0;
        for(int j = 1; j <= r; ++j) 
        { GET(n,j) = t.m(j); }
        for(int j = r+1; j <= NMAX; ++j) 
        { n[j] = 1; }
        init(1);
	}

    Counter& operator++()
	{
        ++ind;
        ++i[1];
        if(i[1] > n[1])
        for(int j = 2; j <= r; ++j)
        {
            i[j-1] = 1;
            ++i[j];
            if(i[j] <= n[j]) break;
        }
        if(i[r] > n[r]) init(0); //set 'done' condition
        return *this;
	}

    bool operator!=(const Counter& other) const
	{
        for(int j = 1; j <= NMAX; ++j)
        { if(i[j] != other.i[j]) return true; }
        return false;
	}
    bool operator==(const Counter& other) const
	{ return !(*this != other); }

    static const Counter *pend;
    static const Counter& done;

    friend inline ostream& operator<<(ostream& s, const Counter& c)
    {
        s << "("; for(int i = 1; i < c.r; ++i){s << c.i[i] << " ";} s << c.i[c.r] << ")";
        return s;
    }
};
#ifdef THIS_IS_MAIN
const Counter* Counter::pend = new Counter;
const Counter& Counter::done(*pend);
#endif                                  

/*
Combine several indices into one, use * to convert tensors efficiently
   \
    \
  ---C====
    /
   /

*/
class Combiner
{
    array<Index,NMAX+1> _leftn; // max dim is 8
    vector<Index> _left1;
    Index _right;
    int _rln; //Number of m>1 'left' indices (indices to be combined into one)
public:
    static Index spec; //Special placeholder index

    //Accessor Methods ----------------------------------------------

    Index right() const { return _right; }
    int rln() const { return _rln; }
    const Index& leftn(int j) const { return _leftn[j]; }

    typedef array<Index,NMAX+1>::const_iterator leftn_it;
    const pair<leftn_it,leftn_it> leftn() const { return make_pair(_leftn.begin()+1,_leftn.begin()+_rln+1); }
    const vector<Index>& left1() const { return _left1; }

    //Constructors --------------------------------------------------

    Combiner() : _rln(0) {}
    Combiner(Index& r, 
	    Index l1 = spec, Index l2 = spec, Index l3 = spec, Index l4 = spec, 
	    Index l5 = spec, Index l6 = spec, Index l7 = spec, Index l8 = spec )
	{
        _leftn[0] = spec;

        //Split given left indices into m==1 and m>1
        _rln = 0;
        if(l1 != spec) { if(l1.m() == 1) _left1.push_back(l1); else _leftn[++_rln] = l1; }
        if(l2 != spec) { if(l2.m() == 1) _left1.push_back(l2); else _leftn[++_rln] = l2; }
        if(l3 != spec) { if(l3.m() == 1) _left1.push_back(l3); else _leftn[++_rln] = l3; }
        if(l4 != spec) { if(l4.m() == 1) _left1.push_back(l4); else _leftn[++_rln] = l4; }
        if(l5 != spec) { if(l5.m() == 1) _left1.push_back(l5); else _leftn[++_rln] = l5; }
        if(l6 != spec) { if(l6.m() == 1) _left1.push_back(l6); else _leftn[++_rln] = l6; }
        if(l7 != spec) { if(l7.m() == 1) _left1.push_back(l7); else _leftn[++_rln] = l7; }
        if(l8 != spec) { if(l8.m() == 1) _left1.push_back(l8); else _leftn[++_rln] = l8; }

        //Set up right index
        int m = 1; foreach(const Index& L, leftn()) m *= L.m();
        r = Index(r.name(),m,r.type(),r.primelevel);
        _right = r;
	}
    
    //Operators -----------------------------------------------------

    friend ITensor operator*(const ITensor& t, const Combiner& c);
    friend inline ITensor operator*(const Combiner& c, const ITensor& t) { return t * c; }

    //Index Methods -------------------------------------------------

    void addleft(Index& r, Index l) 	// Include another left index
	{ 
        if(l.m() == 1) { _left1.push_back(l); return; } else _leftn[++_rln] = l; 
        int m = 1; foreach(const Index& L, leftn()) m *= L.m();
        r = Index(r.name(),m,r.type(),r.primelevel);
        _right = r;
	}
    int findindexn(Index i) const
	{
        for(int j = 1; j <= _rln; ++j)
            if(_leftn[j] == i) return j;
        return 0;
	}
    bool hasindex(Index i) const
	{
        if(i.m() == 1)
        {
            foreach(const Index& L, _left1)
            if(i == L) return true;
            return false;
        }
        for(int j = 1; j <= _rln; ++j) if(_leftn[j] == i) return true;
        return false;
	}

    //Other Methods -------------------------------------------------

    void toITensor(ITensor& res);

    friend inline ostream & operator << (ostream & s, const Combiner & c)
    {
        s << "\nRight index: " << c.right() << "\n";
        s << "Left indices:\n";
        foreach(const Index& l, c.leftn()) s << "	" << l << "\n";
        foreach(const Index& l, c.left1()) s << "	" << l << "\n";
        return s;
    }

}; //class Combiner


enum SweepScheme {ramp_m, fixed_m, fixed_cutoff};

inline void sweepnext(int &l, int &ha, int N, int min_l = 1)
{
    if(ha == 1)
	{
        if(++l == N) 
            l = N-1, ha = 2;
        return;
	}
    if(l-- == min_l) ha = 3;
}

class Sweeps
{
public:
    SweepScheme scheme;
    int Minm;
    vector<int>  Maxm, Niter;
    vector<Real> Cutoff;
    int Nsweep;
    int num_site_center;		// May not be implemented in some cases
    Sweeps(SweepScheme sch, int nsw, int _minm, int _maxm, Real _cut)
	    : scheme(sch), Minm(_minm), Maxm(nsw+1), Niter(nsw+1,4), Cutoff(nsw+1), Nsweep(nsw), num_site_center(2)
	{
        if(scheme == ramp_m)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = _cut; Maxm.at(s) = (int)(_minm + (s-1.0)/nsw * (_maxm - _minm)); }
        }
        else if(scheme == fixed_m || scheme == fixed_cutoff)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = _cut; Maxm.at(s) = _maxm; }
        }

        for(int s = 1; s <= min(Nsweep,4); s++)
        { Niter.at(s) = 10 - s; }
	}
    Real cutoff(int sw) const { return Cutoff.at(sw); }
    int minm(int sw) const { return Minm; }
    int maxm(int sw) const { return Maxm.at(sw); }
    int nsweep() const { return Nsweep; }
    int niter(int sw) const { return Niter.at(sw); }
};

#ifdef COLLECT_PRODSTATS
#define NTIMERS 70
class Prodstats
{
    vector<Real> time;
    vector<int>  tcount;
    vector<cpu_time> cpu;
    vector<bool> timer_running;
public:
    typedef pair<pair<int,int>,int> gitertype;
    map<std::pair<int,int>,int> global;
    map<std::pair<int,int>,int> ps32;
    vector<int> perms_of_3;
    vector<int> perms_of_4;
    vector<int> perms_of_5;
    vector<int> perms_of_6;
    int total, did_matrix;
    int c1,c2,c3,c4;

    Prodstats()
    {
        //for(int i = 0; i <= 20; ++i)
        //for(int j = i; j <= 20; ++j)
            //global[make_pair(j,i)] = 0;
        total = 0;
        did_matrix = 0;
        c1 = c2 = c3 = c4 = 0;
        perms_of_3 = vector<int>(81,0);
        perms_of_4 = vector<int>(256,0);
        perms_of_5 = vector<int>(3125,0);
        perms_of_6 = vector<int>(46656,0);

        time = vector<Real>(NTIMERS,0);
        tcount = vector<int>(NTIMERS,0);
        timer_running = vector<bool>(NTIMERS,false);
        cpu = vector<cpu_time>(NTIMERS);
    }

    void start_section(int j) 
    { 
        if(timer_running.at(j)) Error("Timer already running.");
        timer_running.at(j) = true;
        cpu.at(j) = cpu_time();
    }

    void finish_section(int j)
    {
        if(!timer_running.at(j)) Error("No timer running.");
        timer_running.at(j) = false; 
        cpu_time since(cpu.at(j).sincemark());
        time.at(j) += since.time; 
        tcount.at(j) += 1; 
    }

    void print() const
    {
        cerr << "\n-------- Product Statistics ----------\n";
        cerr << "Global Count: " << endl;
        foreach(gitertype pp, global) cerr << format("(%d,%d) = %d\n")%pp.first.first%pp.first.second%pp.second;
        cerr << "Total = " << total << endl;
        cerr << format("# Matrices = %d (%.2f%%)\n") % did_matrix % (100.0*(1.*did_matrix/(2*total)));

        cerr << "# Case 1 = " << c1 << endl;
        cerr << "# Case 2 = " << c2 << endl;
        cerr << "# Case 3 = " << c3 << endl;
        cerr << "# Case 4 = " << c4 << endl;

        cerr << "Permutations of 3 Count: " << endl;
        for(int j = 0; j < (int) perms_of_3.size(); ++j)
        {
            if(perms_of_3[j] == 0) continue;
            int c = j;
            int i3 = (c%3 == 0 ? 3 : c%3);
            c = (c-i3)/3+1;
            int i2 = (c%3 == 0 ? 3 : c%3);
            c = (c-i2)/3+1;
            int i1 = (c%3 == 0 ? 3 : c%3);
            int idx = ((i1-1)*3+i2-1)*3+i3;
            if(idx != j) cerr << "Incorrect idx val (perms of 3)." << endl;
            cerr << format("(%02d) %d, %d, %d = %d\n") % j % i1 % i2 % i3 % perms_of_3[j];
        }

        cerr << "Permutations of 4 Count: " << endl;
        for(int j = 0; j < (int) perms_of_4.size(); ++j)
        {
            if(perms_of_4[j] == 0) continue;
            int c = j;
            int i4 = (c%4 == 0 ? 4 : c%4);
            c = (c-i4)/4+1;
            int i3 = (c%4 == 0 ? 4 : c%4);
            c = (c-i3)/4+1;
            int i2 = (c%4 == 0 ? 4 : c%4);
            c = (c-i2)/4+1;
            int i1 = (c%4 == 0 ? 4 : c%4);
            int idx = (((i1-1)*4+i2-1)*4+i3-1)*4+i4;
            if(idx != j) cerr << "Incorrect idx val (perms of 4)." << endl;
            cerr << format("(%02d) %d, %d, %d, %d = %d\n") % j % i1 % i2 % i3 % i4 % perms_of_4[j];
        }

        cerr << "Permutations of 5 Count: " << endl;
        for(int j = 0; j < (int) perms_of_5.size(); ++j)
        {
            if(perms_of_5[j] == 0) continue;
            int c = j;
            int i5 = (c%5 == 0 ? 5 : c%5);
            c = (c-i5)/5+1;
            int i4 = (c%5 == 0 ? 5 : c%5);
            c = (c-i4)/5+1;
            int i3 = (c%5 == 0 ? 5 : c%5);
            c = (c-i3)/5+1;
            int i2 = (c%5 == 0 ? 5 : c%5);
            c = (c-i2)/5+1;
            int i1 = (c%5 == 0 ? 5 : c%5);
            int idx = ((((i1-1)*5+i2-1)*5+i3-1)*5+i4-1)*5+i5;
            if(idx != j) cerr << "Incorrect idx val (perms of 5)." << endl;
            cerr << format("(%02d) %d, %d, %d, %d, %d = %d\n") % j % i1 % i2 % i3 % i4 % i5 % perms_of_5[j];
        }

        cerr << "Permutations of 6 Count: " << endl;
        for(int j = 0; j < (int) perms_of_6.size(); ++j)
        {
            //cerr << format("po6[%d] = %d\n") % j % perms_of_6[j];
            if(perms_of_6[j] == 0) continue;
            int c = j;
            int i6 = (c%6 == 0 ? 6 : c%6);
            c = (c-i6)/6+1;
            int i5 = (c%6 == 0 ? 6 : c%6);
            c = (c-i5)/6+1;
            int i4 = (c%6 == 0 ? 6 : c%6);
            c = (c-i4)/6+1;
            int i3 = (c%6 == 0 ? 6 : c%6);
            c = (c-i3)/6+1;
            int i2 = (c%6 == 0 ? 6 : c%6);
            c = (c-i2)/6+1;
            int i1 = (c%6 == 0 ? 6 : c%6);
            int idx = (((((i1-1)*6+i2-1)*6+i3-1)*6+i4-1)*6+i5-1)*6+i6;
            if(idx != j) cerr << "Incorrect idx val (perms of 6)." << endl;
            cerr << format("(%02d) %d, %d, %d, %d, %d, %d = %d\n") % j % i1 % i2 % i3 % i4 % i5 % i6 % perms_of_6[j];
        }

        for(int j = 0; j < (int) time.size(); ++j)
        {
            Real count = tcount.at(j);
            if(time.at(j) > 0) cerr << format("Section %d, Average CPU Time = %.2E\n") % j % (time.at(j)/count);
        }

        for(int j = 0; j < (int) time.size(); ++j)
        {
            if(time.at(j) > 0) cerr << format("Section %d, Total CPU Time = %f\n") % j % time.at(j);
        }
    }
};
#endif
#ifdef COLLECT_PRODSTATS
#define DO_IF_PS(x) { x }
#else
#define DO_IF_PS(x) { }
#endif
#ifdef COLLECT_PRODSTATS
#define START_TIMER(x) { prodstats.start_section(x); }
#else
#define START_TIMER(x) { }
#endif
#ifdef COLLECT_PRODSTATS
#define STOP_TIMER(x) { prodstats.finish_section(x); }
#else
#define STOP_TIMER(x) { }
#endif

class Measure
{
public:
    vector<Real> dat;
    Measure() { dat.reserve(100); }

    void putin(Real x) { dat.push_back(x); }

    Real ave() const
	{
        if(dat.empty()) return 0;
        Real s = 0;
        foreach(Real d, dat) s += d;
        return s/dat.size();
	}

    Real sigma() const
	{
        if(dat.size() < 2) return 0;
        Real av = ave(), s = 0;
        foreach(Real d, dat) s += sqr(d - av);
        return sqrt(s/dat.size());
	}
    Real err() const
	{
        Real n = dat.size();
        if(n > 1) return sigma()/sqrt(n-1);
        else return sigma()/sqrt(n);
	}
    Real binerr(int bs) const
	{
        int n = dat.size();
        if(n < bs) return 0;
        int nbin = n / bs;
        Vector bin(nbin);
        bin = 0.0;
        for(int i = 0; i < n; i++)
        {
            if(i/bs > (nbin-1)) break;
            bin.el(i/bs) += GET(dat,i);
        }
        bin *= 1.0/bs;
        Measure bmeas;
        for(int i = 1; i <= nbin; ++i) bmeas.putin(bin(i));
        return bmeas.err();
	}
    Real bootstrap(int nresamples = 1000)
    { //Uses the bootstrap procedure to estimate the std deviation
        Measure boot;
        int n = dat.size();
        if(n < 2) return 0.0;
        for(int sample = 1; sample <= nresamples; ++sample)
        {
            Real avg = 0.0;
            for(int i = 1; i <= n; ++i)
            {
                int which = int(ran1()*n);
                avg += dat.at(which);
            }
            avg /= n;
            boot.putin(avg);
        }
        return boot.sigma();
    }

    Real bootstrap_2c(Real prefactor,Measure firstc, bool correlated = true, int nresamples = 1000)
    { //Uses the bootstrap procedure to estimate the second cumulant error
        Measure boot2c;
        Real n = dat.size();
        if(n < 2) return 0.0;
        if(firstc.dat.size() != n) error("Measurements don't have the same size data sets.");
        for(int sample = 1; sample <= nresamples; ++sample)
        {
            Real avg2 = 0.0;
            Real avg = 0.0;
            int which;
            for(int i = 1; i <= n; ++i)
            { //Sample with replacement (replacement means 'which' can take the same value more than once)
                which = int(ran1()*n);
                avg2 += dat[which];
                if(!correlated) which = int(ran1()*n); //use <H^2> and <H> from the same METTS
                avg += firstc.dat[which];
            }
            avg2 /= n; avg /= n;
            boot2c.putin(prefactor*(avg2-avg*avg));
        }
        return boot2c.sigma();
    }

};


inline void System(const char* cstr) { int res = system(cstr); res++; }
inline void System(const string str) { System(str.c_str()); }
inline void System(const stringstream& s) { System(s.str()); }
inline void System(const format fmt) { System(fmt.str()); }

//Vector data types; x, y and error data
inline void writedata(const char* cstr, const Vector& xdat,const Vector& ydat, const Vector& edat, bool do_plot_self = false)
{
    ofstream f(cstr);
    if(xdat.Length() != ydat.Length()) Error("xdat and ydat Lengths don't match.");
    for(int j = 1; j <= xdat.Length(); ++j) f << format("%.10f %.15f %.15f\n") % xdat(j) % ydat(j) % edat(j);
    f.close();
    //if(do_plot_self) System(format("$HOME/tools/plot_self %s") % cstr);
    if(do_plot_self) System(format("plot_self %s") % cstr);
}
inline void writedata(const string str, const Vector& xdat, const Vector& ydat, const Vector& edat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,edat,do_plot_self); }
inline void writedata(const stringstream& s,const Vector& xdat, const Vector& ydat, const Vector& edat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,edat,do_plot_self); }
inline void writedata(const format fmt, const Vector& xdat, const Vector& ydat, const Vector& edat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,edat,do_plot_self); }

//Vector data types; x, y data
inline void writedata(const char* cstr, const Vector& xdat,const Vector& ydat, bool do_plot_self = false)
{
    Vector edat = ydat; edat = 0;
    writedata(cstr,xdat,ydat,edat,do_plot_self);
}
inline void writedata(const string str, const Vector& xdat, const Vector& ydat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,do_plot_self); }
inline void writedata(const stringstream& s,const Vector& xdat, const Vector& ydat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,do_plot_self); }
inline void writedata(const format fmt, const Vector& xdat, const Vector& ydat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,do_plot_self); }

//std vector data types; x, y and error data
template<typename T>
void writedata(const char* cstr, const vector<T>& xdat,const vector<T>& ydat, const vector<T>& edat, bool do_plot_self = false)
{
    Vector X((int) xdat.size()); for(int n = 1; n <= X.Length(); ++n) X(n) = xdat.at(n-1);
    Vector Y((int) ydat.size()); for(int n = 1; n <= Y.Length(); ++n) Y(n) = ydat.at(n-1);
    Vector E((int) edat.size()); for(int n = 1; n <= E.Length(); ++n) E(n) = edat.at(n-1);
    writedata(cstr,X,Y,E,do_plot_self);
}
template<typename T>
inline void writedata(const string str, const vector<T>& xdat, const vector<T>& ydat, const vector<T>& edat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,edat,do_plot_self); }
template<typename T>
inline void writedata(const stringstream& s,const vector<T>& xdat, const vector<T>& ydat, const vector<T>& edat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,edat,do_plot_self); }
template<typename T>
inline void writedata(const format fmt, const vector<T>& xdat, const vector<T>& ydat, const vector<T>& edat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,edat,do_plot_self); }

//std vector data types; x, y and error data
template<typename T>
void writedata(const char* cstr, const vector<T>& xdat,const vector<T>& ydat, bool do_plot_self = false)
{
    Vector X((int) xdat.size()); for(int n = 1; n <= X.Length(); ++n) X(n) = xdat.at(n-1);
    Vector Y((int) ydat.size()); for(int n = 1; n <= Y.Length(); ++n) Y(n) = ydat.at(n-1);
    writedata(cstr,X,Y,do_plot_self);
}
template<typename T>
inline void writedata(const string str, const vector<T>& xdat, const vector<T>& ydat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,do_plot_self); }
template<typename T>
inline void writedata(const stringstream& s,const vector<T>& xdat, const vector<T>& ydat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,do_plot_self); }
template<typename T>
inline void writedata(const format fmt, const vector<T>& xdat, const vector<T>& ydat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,do_plot_self); }

template<typename T>
void writedata(const char* cstr,const vector<T>& ydat, bool do_plot_self = false)
{
    Vector X((int) ydat.size()); for(int n = 1; n <= X.Length(); ++n) X(n) = n-1;
    Vector Y((int) ydat.size()); for(int n = 1; n <= Y.Length(); ++n) Y(n) = ydat.at(n-1);
    writedata(cstr,X,Y,do_plot_self);
}
template<typename T>
inline void writedata(const string str, const vector<T>& ydat, bool do_plot_self = false) { writedata(str.c_str(),ydat,do_plot_self); }
template<typename T>
inline void writedata(const stringstream& s, const vector<T>& ydat, bool do_plot_self=false) { writedata(s.str(),ydat,do_plot_self); }
template<typename T>
inline void writedata(const format fmt, const vector<T>& ydat, bool do_plot_self=false) { writedata(fmt.str(),ydat,do_plot_self); }

inline void writedata(const char* cstr, const Vector& dat, Real Delta = 1, bool do_plot_self = false)
{
    ofstream f(cstr);
    for(int j = 1; j <= dat.Length(); ++j) f << format("%.10f %.10f\n") % (Delta*j) % dat(j);
    f.close();
    //if(do_plot_self) System(format("$HOME/tools/plot_self %s") % cstr);
    if(do_plot_self) System(format("plot_self %s") % cstr);
}
inline void writedata(const string str, const Vector& dat, Real Delta = 1, bool do_plot_self = false) { writedata(str.c_str(),dat,Delta,do_plot_self); }
inline void writedata(const stringstream& s,const Vector& dat, Real Delta = 1, bool do_plot_self=false) { writedata(s.str(),dat,Delta,do_plot_self); }
inline void writedata(const format fmt, const Vector& dat, Real Delta = 1, bool do_plot_self=false) { writedata(fmt.str(),dat,Delta,do_plot_self); }

inline void writedata(const char* cstr, const Matrix& dat, bool do_plot_self = false)
{           
    ofstream f(cstr);
    for(int r = 1; r <= dat.Nrows(); ++r)
    for(int c = 1; c <= dat.Ncols(); ++c) 
        f << r SP c SP dat(r,c) << "\n";
    f.close();
    //if(do_plot_self) System(format("$HOME/tools/plot_self %s") % cstr);
    if(do_plot_self) System(format("plot_self %s") % cstr);
}           
inline void writedata(const string str, const Matrix& dat, bool do_plot_self=false) { writedata(str.c_str(),dat,do_plot_self); }
inline void writedata(const stringstream& s,const Matrix& dat, bool do_plot_self=false) { writedata(s.str(),dat,do_plot_self); }
inline void writedata(const format fmt, const Matrix& dat, bool do_plot_self=false) { writedata(fmt.str(),dat,do_plot_self); }




#ifdef THIS_IS_MAIN
void reportnew() {}
//Real ran1();
Real ran1(int);

//int Internal::IndexDat::indcount = 0;
UniqueID Internal::IndexDat::lastID; 
bool printdat = false;
bool writeops = false;
ofstream big_op_file;
Internal::IndexDat IndexDatNull(makeNull);
Internal::IndexDat IndReDat(makeReIm);
Internal::IndexDat IndReDatP(makeReImP);
Internal::IndexDat IndReDatPP(makeReImPP);
Internal::IndexDat IndEmptyVDat(makeEmptyV);
Index Combiner::spec = Index("spec",2);
Index IndNull(makeNull);
Index IndReIm(makeReIm);
Index IndReImP(makeReImP);
Index IndReImPP(makeReImPP);
Index IndEmptyV(makeEmptyV);
IndexVal IVNull(IndNull,1);
ITensor Complex_1(makeComplex_1), Complex_i(makeComplex_i), ConjTensor(makeConjTensor);
Vector lastd(1);
int newtotalsize = 0;

//Debugging and profiling stuff:
bool catch_debug = false;
#ifdef COLLECT_PRODSTATS
Prodstats prodstats;
#endif
#ifdef COUNT_COPIES
int copycount = 0;
#endif

#endif //end ifdef THIS_IS_MAIN

#endif
