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



///
/// Index: represents a tensor index of fixed bond dimension m.
/// The == operator can be used to determine if two Index's match
/// (copies of the same Index instance).
/// An Index can be made temporarily distinct from other copies 
/// by increasing its "primelevel".
///
class Index
    {
    public:

    /// Returns the bond dimension of an Index.
    /// @return Bond dimension of this Index
    int 
    m() const;

    const boost::uuids::uuid&
    Ind() const;

    /// Returns the IndexType of this Index.
    /// @return IndexType of this Index
    IndexType 
    type() const;

    /// Returns the name of this Index.
    /// @return Name of this Index
    std::string 
    name() const;

    /// Returns the name of this Index with primes removed.
    /// @return Name of this Index without primes
    const std::string&
    rawname() const;

    /// Sets the name of this Index.
    void 
    setname(const std::string& newname);

    /// Returns a string version of this Index's bond dimension.
    /// @return String version of bond dimension
    std::string 
    showm() const;

    /// Returns a unique Real number identifying this Index.
    /// Useful for rapidly checking that two Index instances match.
    /// @return Unique Real id number
    Real 
    uniqueReal() const;

    /// Returns true if Index was default initialized.
    /// @return True if default initialized
    bool 
    isNull() const;
    /// Returns true if Index was NOT default initialized.
    /// @return True if not default initialized
    bool 
    isNotNull() const;

    /// Returns the number of copies of this Index.
    /// @return Number of copies of Index
    int 
    count() const;

    /// Returns the prime level of this Index
    /// @return Prime level
    int 
    primeLevel() const;
    /// Sets the prime level to a specified value.
    /// @param plev New prime level
    void 
    primeLevel(int plev);

    /// Returns the Arrow direction of this Index.
    /// @return Arrow direction
    Arrow 
    dir() const { return Out; }

    //-----------------------------------------------
    //Index: Constructors

    /// Default constructor.
    /// For default constructed Index's, isNull() returns true.
    Index();

    /// Index constructor.
    /// @param name Name of Index for printing purposes
    /// @param mm   Bond dimension 
    /// @param it   IndexType of Index (Link, Site)
    /// @param plev Prime level
    Index(const std::string& name, 
          int mm = 1, 
          IndexType it=Link, 
          int plev = 0);

    /// Input stream constructor.
    /// @param s Binary input stream
    Index(std::istream& s) { read(s); }

    /// Prime level copy constructor.
    /// @param pt       PrimeType to use
    /// @param other    Index to copy and prime
    /// @param primeinc Prime increment level
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
