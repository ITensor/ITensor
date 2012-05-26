//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEX_H
#define __ITENSOR_INDEX_H
#include <string>
#include "types.h"
#include "boost/intrusive_ptr.hpp"
#include "boost/uuid/uuid.hpp"
#include "boost/uuid/random_generator.hpp"
#include "boost/uuid/string_generator.hpp"

enum Arrow { In = -1, Out = 1 };

Arrow inline
operator*(const Arrow& a, const Arrow& b)
    { return (int(a)*int(b) == In) ? In : Out; }

const Arrow Switch = In*Out;

inline std::ostream& 
operator<<(std::ostream& s, Arrow D)
    { if(D == In) s << "In"; else s << "Out"; return s; }

enum IndexType { Link, Site, ReIm };
static const char* indextypename[] = { "Link","Site","ReIm" };

enum PrimeType { primeLink, primeSite, primeBoth, primeNone };

inline std::ostream& 
operator<<(std::ostream& s, const IndexType& it)
    { 
    if(it == Link) s << "Link"; 
    else if(it == Site) s << "Site"; 
    else if(it == ReIm) s << "ReIm"; 
    return s; 
    }

inline int 
IndexTypeToInt(IndexType it)
    {
    if(it == Link) return 1;
    if(it == Site) return 2;
    if(it == ReIm) return 3;
    Error("No integer value defined for IndexType.");
    return -1;
    }
inline IndexType 
IntToIndexType(int i)
    {
    if(i == 1) return Link;
    if(i == 2) return Site;
    if(i == 3) return ReIm;
    std::cout << boost::format("No IndexType value defined for i=%d\n")%i 
              << std::endl;
    Error("Undefined IntToIndexType value");
    return Link;
    }

std::string inline
putprimes(std::string s, int plev = 0)
    { 
    for(int i = 1; i <= plev; ++i) 
        s += "\'"; 
    return s;
    }

std::string inline
nameindex(IndexType it, int plev = 0)
    { 
    return putprimes(std::string(indextypename[(int)it]),plev); 
    }

std::string inline
nameint(std::string f,int ix)
    { 
    std::stringstream ss; 
    ss << f << ix; 
    return ss.str(); 
    }

enum Imaker {makeReIm,makeReImP,makeReImPP,makeNull};

#define UID_NUM_PRINT 2
inline std::ostream& 
operator<<(std::ostream& s, const boost::uuids::uuid& id)
    { 
    s.width(2);
    for(boost::uuids::uuid::size_type i = id.size()-UID_NUM_PRINT; i < id.size(); ++i) 
        {
        s << static_cast<unsigned int>(id.data[i]);
        }
    s.width(0);
    return s; 
    }

struct UniqueID
    {
    boost::uuids::uuid id;

    UniqueID() : id(boost::uuids::random_generator()()) { }

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

    operator boost::uuids::uuid() const { return id; }

    friend inline std::ostream&
    operator<<(std::ostream& s, const UniqueID& uid) 
        { 
        s << uid.id; 
        return s; 
        }
    };

class IndexDat;
struct IndexVal;



//
// Index
//

class Index
    {
    public:

    int 
    m() const;

    const boost::uuids::uuid&
    Ind() const;

    IndexType 
    type() const;

    std::string 
    name() const;

    const std::string&
    rawname() const;

    void 
    setname(const std::string& newname);

    std::string 
    showm() const;

    Real 
    uniqueReal() const;

    bool 
    isNull() const;
    bool 
    isNotNull() const;

    int 
    count() const;

    int 
    primeLevel() const;
    void 
    primeLevel(int plev);

    Arrow 
    dir() const { return Out; }

    //-----------------------------------------------
    //Index: Constructors

    Index();

    Index(const std::string& name, int mm = 1, IndexType it=Link, int plev = 0);

    Index(std::istream& s) { read(s); }

    Index(PrimeType pt,const Index& other, int primeinc = 1);

    static const Index& Null()
        {
        static const Index Null_(makeNull);
        return Null_;
        }

    static const Index& IndReIm()
        {
        static const Index IndReIm_(makeReIm);
        return IndReIm_;
        }

    static const Index& IndReImP()
        {
        static const Index IndReImP_(makeReImP);
        return IndReImP_;
        }

    static const Index& IndReImPP()
        {
        static const Index IndReImPP_(makeReImPP);
        return IndReImPP_;
        }

    //-----------------------------------------------
    //Index: Operators

    // rel_ops defines the other comparisons based on == and <
    bool 
    operator==(const Index& other) const;

    bool 
    noprime_equals(const Index& other) const;

    bool 
    operator<(const Index& other) const;

    IndexVal operator()(int i) const;


    //-----------------------------------------------
    //Index: Prime methods

    void 
    mapprime(int plevold, int plevnew, PrimeType pr = primeBoth);

    void 
    doprime(PrimeType pr, int inc = 1);

    Index 
    primed(int inc = 1) const { return Index(primeBoth,*this,inc); }

    friend inline Index
    primed(const Index& I, int inc = 1) { return Index(primeBoth,I,inc); }

    Index 
    deprimed() const;

    void 
    noprime(PrimeType pt = primeBoth) { doprime(pt,-primelevel_); }

    friend std::ostream& 
    operator<<(std::ostream & s, const Index & t);

    //-----------------------------------------------
    //Index: Other methods

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    void 
    print(std::string name = "") const
        { std::cerr << "\n" << name << " =\n" << *this << "\n"; }

    void 
    conj() { } //for forward compatibility with arrows

    private:

    friend class IQIndex;

    //void 
    //set_m(int newm);

    //Constructor for static Index's
    explicit
    Index(Imaker im);

    /////////////
    //
    // Data Members

    boost::intrusive_ptr<IndexDat> p;

    int primelevel_; 

    //
    /////////////

    }; //class Index


//
// IndexDat
//

//Storage for Indexes
class IndexDat
    {
    public:


    //////////////
    //
    // Public Data Members

    IndexType _type;
    boost::uuids::uuid ind;
    int m_;
    Real ur;
    std::string sname;

    //
    //////////////

    void 
    setUniqueReal();

    IndexDat(const std::string& name="", int mm = 1,IndexType it=Link);

    //For use with read/write functionality of Index class
    IndexDat(const std::string& ss, int mm, IndexType it, const boost::uuids::uuid& ind_);


    static IndexDat* 
    Null()
        {
        static IndexDat Null_(makeNull);
        return &Null_;
        }

    static IndexDat* 
    ReImDat()
        {
        static IndexDat ReImDat_(makeReIm);
        return &ReImDat_;
        }

    static IndexDat* 
    ReImDatP()
        {
        static IndexDat ReImDatP_(makeReImP);
        return &ReImDatP_;
        }

    static IndexDat* 
    ReImDatPP()
        {
        static IndexDat ReImDatPP_(makeReImPP);
        return &ReImDatPP_;
        }

    friend void 
    intrusive_ptr_add_ref(IndexDat* p);

    friend void 
    intrusive_ptr_release(IndexDat* p);

    int 
    count() const { return numref; }

    private:

    //////////////
    //
    // (Private) Data Members
    
    mutable unsigned int numref;

    const bool is_static_;

    //
    /////////////

    explicit
    IndexDat(Imaker im);

    static const UniqueID& 
    nextID();

    //These constructors are not implemented
    //to disallow copying
    IndexDat(const IndexDat&);
    void operator=(const IndexDat&);

    }; //class IndexDat


//
// IndexVal
//

struct IndexVal
    {
    Index ind; 
    int i;

    IndexVal();

    IndexVal(const Index& index, int i_);

    static const IndexVal& 
    Null()
        {
        static const IndexVal Null_(makeNull);
        return Null_;
        }

    bool 
    operator==(const IndexVal& other) const; 

    friend IndexVal
    primed(const IndexVal& iv, int inc = 1);

    friend std::ostream& 
    operator<<(std::ostream& s, const IndexVal& iv);

    private:

    explicit
    IndexVal(Imaker im);

    };

template <class T> T 
conj(T res) 
    { 
    res.conj(); 
    return res; 
    }

#endif
