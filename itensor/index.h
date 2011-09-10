#ifndef __ITENSOR_INDEX_H
#define __ITENSOR_INDEX_H
#include <string>
#include "types.h"
#include "boost/array.hpp"
#include "boost/format.hpp"
#include "boost/intrusive_ptr.hpp"
#include "boost/uuid/uuid.hpp"
#include "boost/uuid/random_generator.hpp"
#include "boost/uuid/string_generator.hpp"
using boost::intrusive_ptr;
using boost::uuids::uuid;
using boost::uuids::random_generator;
using boost::uuids::string_generator;

enum Arrow { In = -1, Out = 1 };

inline Arrow operator*(const Arrow& a, const Arrow& b)
{ return (int(a)*int(b) == In) ? In : Out; }

const Arrow Switch = In*Out;

inline std::ostream& operator<<(std::ostream& s, const Arrow& D)
{ if(D == In) s << "In"; else s << "Out"; return s; }

enum IndexType { Link, Site, ReIm, Virtual };
static const char * indextypename[] = { "Link","Site","ReIm","Virtual" };

enum PrimeType { primeLink, primeSite, primeBoth, primeNone };

inline std::ostream& operator<<(std::ostream& s, const IndexType& it)
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
    cerr << boost::format("No IndexType value defined for i=%d\n")%i,Error("");
    return Virtual;
}

inline std::string putprimes(std::string s, int plev = 0)
{ for(int i = 1; i <= plev; ++i) s += "\'"; return s;}

inline std::string nameindex(IndexType it, int plev = 0)
{ return putprimes(std::string(indextypename[(int)it]),plev); }

inline std::string nameint(std::string f,int ix)
{ std::stringstream ss; ss << f << ix; return ss.str(); }

enum Imaker {makeReIm,makeReImP,makeReImPP,makeEmptyV,makeNull};

#define UID_NUM_PRINT 2
inline std::ostream& operator<<(std::ostream& s, const uuid& id)
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

    friend inline std::ostream& operator<<(std::ostream& s, const UniqueID& uid) { s << uid.id; return s; }
};

namespace {
inline int prime_number(int n)
{
    static const boost::array<int,54> plist = { { 
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 
    79, 83, 89, 97, 101, 103, 107, 109, 113, 
    127, 131, 137, 139, 149, 151, 157, 163, 
    167, 173, 179, 181, 191, 193, 197, 199, 
    211, 223, 227, 229, 233, 239, 241, 251 
    } };
    return plist.at(n);
}
}
//Storage for Index's
class IndexDat
{
    mutable unsigned int numref;
    const bool is_static_;
public:
    static UniqueID lastID;

    IndexType _type;
    uuid ind;
    int m_;
    Real ur;
    std::string sname;

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
	}

    IndexDat(std::string name="", int mm = 1,IndexType it=Link) :
    numref(0), is_static_(false),
    _type(it), 
    ind(++lastID),
    m_(mm), 
    sname(name)
	{ 
        if(it == ReIm) Error("bad call to create IndexDat with type ReIm");
        assert((it==Virtual ? (mm==1) : true)); //If type is Virtual, m must be 1
        set_unique_Real();
	}

    //For use with read/write functionality of Index class
    IndexDat(std::string ss, int mm, IndexType it, uuid ind_) :
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
    m_( (im==makeNull || im==makeEmptyV) ? 1 : 2)
	{ 
        string_generator gen;
        if(im==makeNull || im==makeEmptyV) 
        { ind = gen("{00000000-0000-0000-0000-000000000000}"); }
        else                               
        { ind = gen("{10000000-0000-0000-0000-000000000000}"); }

        if(im == makeNull)
        {
            sname = "Null";
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
private:
    IndexDat(const IndexDat&);
    void operator=(const IndexDat&);
};

extern IndexDat IndexDatNull, IndReDat, IndReDatP, IndReDatPP, IndEmptyVDat;

struct IndexVal;

class Index
{
private:
    intrusive_ptr<IndexDat> p;
protected:
    void set_m(int newm) { p->m_ = newm; }
public:
    int primelevel; 

    inline int m() const { return p->m_; }
    uuid Ind() const { return p->ind; }
    inline IndexType type() const { return p->_type; }

    std::string name() const  { return putprimes(rawname(),primelevel); }
    std::string rawname() const { return p->sname; }
    void setname(std::string newname) { p->sname = newname; }

    inline std::string showm() const { return (boost::format("m=%d")%(p->m_)).str(); }
    Real unique_Real() const { assert(p!=0); return p->ur*(1+primelevel); }
    inline bool is_null() const { return (p == &IndexDatNull); }
    inline bool is_not_null() const { return (p != &IndexDatNull); }
    int count() const { return p->count(); }

    void setPrimeLevel(int plev) { primelevel = plev; }

    //-----------------------------------------------
    //Index: Constructors

    Index() : p(&IndexDatNull), primelevel(0) { }

    Index(std::string name, int mm = 1, IndexType it=Link, int plev = 0) 
	: p(new IndexDat(name,mm,it)), primelevel(plev) { }

    Index(std::istream& s) { read(s); }

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

    //-----------------------------------------------
    //Index: Operators

    // rel_ops defines the other comparisons based on == and <
    bool operator==(const Index& other) const 
	{ return unique_Real() == other.unique_Real(); }

    bool operator<(const Index& other) const 
	{ return (unique_Real() < other.unique_Real()); }

    IndexVal operator()(int i) const;

    bool noprime_equals(const Index& other) const
	{ return (p->ur == other.p->ur); }

    //-----------------------------------------------
    //Index: Prime methods

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

    friend inline std::ostream & operator<<(std::ostream & s, const Index & t)
    {
        if(t.name() != "" && t.name() != " ") s << t.name() << "/";
        return s << nameindex(t.type(),t.primelevel) << "-" << t.Ind() << ":" << t.m();
    }

    //-----------------------------------------------
    //Index: Other methods

    void write(std::ostream& s) const 
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

    void read(std::istream& s)
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
        std::string ss(newname); delete newname;
        p = new IndexDat(ss,mm,IntToIndexType(t),ind);
    }

    void print(std::string name = "") const
    { cerr << "\n" << name << " =\n" << *this << "\n"; }

    void conj() { } //for forward compatibility with arrows

}; //class Index
extern Index IndNull, IndReIm, IndReImP, IndReImPP, IndEmptyV;

template <class T> 
T conj(T res) { res.conj(); return res; }

struct IndexVal
{
    Index ind; 
    int i;
    IndexVal() : ind(IndNull),i(0) { }
    IndexVal(const Index& index, int i_) : ind(index),i(i_) { assert(i <= ind.m()); }
    bool operator==(const IndexVal& other) const { return (ind == other.ind && i == other.i); }
    inline friend std::ostream& operator<<(std::ostream& s, const IndexVal& iv)
    { return s << "IndexVal: i = " << iv.i << ", ind = " << iv.ind << "\n"; }
    IndexVal primed() const { return IndexVal(ind.primed(),i); }
};
extern IndexVal IVNull;

inline IndexVal Index::operator()(int i) const 
{ return IndexVal(*this,i); }

#ifdef THIS_IS_MAIN
UniqueID IndexDat::lastID; 
IndexDat IndexDatNull(makeNull);
IndexDat IndReDat(makeReIm);
IndexDat IndReDatP(makeReImP);
IndexDat IndReDatPP(makeReImPP);
IndexDat IndEmptyVDat(makeEmptyV);
Index IndNull(makeNull);
Index IndReIm(makeReIm);
Index IndReImP(makeReImP);
Index IndReImPP(makeReImPP);
Index IndEmptyV(makeEmptyV);
IndexVal IVNull(IndNull,1);
#endif //THIS_IS_MAIN

#endif
