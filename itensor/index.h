//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEX_H
#define __ITENSOR_INDEX_H
#include "itensor/global.h"
#include "itensor/indextype.h"
#include "itensor/indexname.h"
#include "itensor/arrow.h"

namespace itensor {

//Forward declarations
class IndexVal;

//
// Index
//
// Represents a tensor index of fixed bond dimension m.
//
// Can be compared with == operator (returns true if both
// are copies of the same Index instance).
//
// To make an Index distinct from other copies, increase its primeLevel.
//
class Index
    {
    public:
    using IDGenerator = std::mt19937;
    using id_type = IDGenerator::result_type;
    using indexval_type = IndexVal;
    using prime_type = int;
    using extent_type = int;
    private:
    id_type id_;
    prime_type primelevel_; 
    extent_type m_;
    IndexType type_;
    IndexName name_;
    public:

    //
    // Constructors
    //

    Index();

    // Name of Index is used for printing purposes
    explicit
    Index(const std::string& name, 
          long m = 1, 
          IndexType it = Link, 
          int primelev = 0);

    //
    // Accessor Methods
    //

    // Returns the bond dimension
    long 
    m() const { return m_; }

    // Returns the prime level
    int 
    primeLevel() const { return primelevel_; }
    // Sets the prime level to a specified value.
    Index& 
    primeLevel(int plev);

    // Returns the IndexType
    IndexType
    type() const { return type_; }

    // Returns the name of this Index
    std::string 
    name() const;

    // Returns the name of this Index with primes removed
    std::string
    rawname() const { return std::string(name_.c_str()); }

    // Evaluates to false if Index is default constructed.
    explicit operator bool() const;

    // (Explicitly) convertible to integer types
    explicit operator int() const { return m(); }
    explicit operator long() const { return m(); }
    explicit operator size_t() const { return m(); }


    // Returns the Arrow direction of this Index
    Arrow 
    dir() const { return Out; }
    void 
    dir(Arrow ndir) const {  }


    //
    // Prime level methods
    //

    // Increase primelevel by 1 (or by optional amount inc)
    Index& 
    prime(int inc = 1);

    // Increase primelevel by 1 (or optional amount inc)
    // if type matches this Index or type==All
    Index& 
    prime(IndexType type, int inc = 1);

    // Set primelevel to zero (optionally only if type matches)
    Index& 
    noprime(IndexType type = All) { prime(type,-primelevel_); return *this; }

    // Switch primelevel from plevold to plevnew
    // Has no effect if plevold doesn't match current primelevel
    Index& 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //
    // Operators
    //

    // Equality (==) operator
    // True if other Index is a copy of this Index and has same primelevel
    bool 
    operator==(const Index& other) const;

    bool 
    operator!=(const Index& other) const { return !operator==(other); }

    // Check if other Index is a copy of this, ignoring primeLevel
    bool 
    noprimeEquals(const Index& other) const;

    // Useful for sorting Index objects
    bool 
    operator>(const Index& other) const;
    bool 
    operator<(const Index& other) const;

    // Creates an IndexVal from this Index with value i
    IndexVal 
    operator()(long i) const;

    // Creates an IndexVal from this Index 
    // with prime level plev and value (nplev-plev)
    // (when passed to prime function, indicates map from plev->nplev)
    IndexVal 
    operator()(long plev, long nplev) const;

    //
    // Other methods
    //

    // Conjugate this Index.
    // Currently has no effect; exists for forward compatibility
    // with Arrows and interface compatibility with class IQIndex.
    void 
    dag() { } //for forward compatibility with arrows

    // Write Index to binary output stream.
    void 
    write(std::ostream& s) const;

    // Read Index from binary input stream.
    Index& 
    read(std::istream& s);

    id_type
    id() const { return id_; }


    private:

    /////////////
    /////////////

    Index::id_type 
    generateID()
        {
        static Index::IDGenerator rng(std::time(NULL) + getpid());
        return rng();
        }

    }; //class Index


//
// IndexVal
//
// Class pairing an Index of dimension m
// with a specific value i where 1 <= i <= m
//
class IndexVal
    {
    public:

    Index index;
    long val;

    IndexVal();

    IndexVal(const Index& index, long val_);

    bool
    operator==(const IndexVal& other) const;

    bool
    operator!=(const IndexVal& other) const { return !operator==(other); }

    bool
    operator==(const Index& I) const { return index == I; }

    long
    m() const { return index.m(); }

    explicit operator bool() const { return bool(index); }

    IndexVal& 
    prime(int inc = 1) { index.prime(inc); return *this; }

    IndexVal& 
    prime(IndexType type, int inc = 1) { index.prime(type,inc); return *this; }

    IndexVal& 
    noprime(IndexType type = All) { index.noprime(type); return *this; }

    IndexVal& 
    mapprime(int plevold, int plevnew, IndexType type = All) 
        { index.mapprime(plevold,plevnew,type); return *this; }

    void
    dag() { }

    };

bool inline
operator==(const Index& I, const IndexVal& iv)
    {
    return iv.operator==(I);
    }

Index inline
dag(Index res) { res.dag(); return res; }

IndexVal inline
dag(IndexVal res) { res.dag(); return res; }

template<typename... VarArgs>
Index
prime(Index I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
Index
noprime(Index I, VarArgs&&... vargs) { I.noprime(std::forward<VarArgs>(vargs)...); return I; }

//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
Index inline
mapprime(Index I, int plevold, int plevnew, IndexType type = All)
    { I.mapprime(plevold,plevnew,type); return I; }

template<typename... VarArgs>
IndexVal
prime(IndexVal I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
IndexVal
noprime(IndexVal I, VarArgs&&... vargs) { I.noprime(std::forward<VarArgs>(vargs)...); return I; }

//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
IndexVal inline
mapprime(IndexVal I, int plevold, int plevnew, IndexType type = All)
    { I.mapprime(plevold,plevnew,type); return I; }


//Returns a string version of this Index's bond dimension.
std::string
showm(const Index& I);

std::string 
nameint(const std::string& f, int n);

std::ostream& 
operator<<(std::ostream & s, const Index &t);

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv);


//
//
// Implementations
//
//

inline Index::
Index() 
    : 
    id_(0),
    primelevel_(0),
    m_(1),
    type_(NullInd)
    { }

inline Index::
Index(const std::string& name, long m, IndexType type, int plev) 
    : 
    id_(generateID()),
    primelevel_(plev),
    m_(m),
    type_(type),
    name_(name.c_str())
    { 
#ifdef DEBUG
    if(type_ == All) Error("Constructing Index with type All disallowed");
    if(type_ == NullInd) Error("Constructing Index with type NullInd disallowed");
#endif
    }

bool inline Index::
operator==(const Index& other) const 
    { 
    return (id_ == other.id_) && (primelevel_ == other.primelevel_); 
    }

bool inline Index::
noprimeEquals(const Index& other) const
    { 
    return (id_ == other.id_);
    }

bool inline Index::
operator>(const Index& other) const 
    { 
    if(m_ == other.m_) 
        {
        if(id_ == other.id_) return primelevel_ > other.primelevel_;
        return id_ > other.id_;
        }
    return m_ > other.m_;
    }

bool inline Index::
operator<(const Index& other) const
    {
    if(m_ == other.m_) 
        {
        if(id_ == other.id_) return primelevel_ < other.primelevel_;
        return id_ < other.id_;
        }
    return m_ < other.m_;
    }

IndexVal inline Index::
operator()(long val) const { return IndexVal(*this,val); }

IndexVal inline Index::
operator()(long plev, long nplev) const { return IndexVal(itensor::prime(*this,plev),nplev - plev); }

inline
Index& Index::
primeLevel(int plev) 
    { 
    primelevel_ = plev; 
#ifdef DEBUG
    if(primelevel_ < 0)
        Error("Negative primeLevel");
#endif
    return *this;
    }

inline
Index& Index::
prime(int inc) 
    { 
    primelevel_ += inc; 
#ifdef DEBUG
    if(primelevel_ < 0)
        {
        Error("Negative primeLevel");
        }
#endif
    return *this;
    }

void
add(Args& args, 
    const Args::Name& name, 
    IndexType it);

IndexType
getIndexType(const Args& args, 
             const Args::Name& name);

IndexType
getIndexType(const Args& args, 
             const Args::Name& name, 
             IndexType default_val);

} //namespace itensor

#endif
