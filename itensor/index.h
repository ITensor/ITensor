//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEX_H
#define __ITENSOR_INDEX_H
#include "global.h"
#include "boost/shared_ptr.hpp"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

enum IndexType { Link, Site, All };

//Forward declarations
struct IndexDat;
class IndexVal;

typedef boost::shared_ptr<IndexDat>
IndexDatPtr;

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
    //
    // Constructors
    //

    // For a default constructed Index, isNull() returns true.
    Index();

    // Name of Index is used for printing purposes
    explicit
    Index(const std::string& name, 
          int m = 1, 
          IndexType it = Link, 
          int primelev = 0);

    //
    // Accessor Methods
    //

    // Returns the bond dimension
    int 
    m() const;

    // Returns the prime level
    int 
    primeLevel() const;
    // Sets the prime level to a specified value.
    Index& 
    primeLevel(int plev);

    // Returns a unique Real number identifying this Index.
    // Useful for efficiently checking that sets of indices match.
    Real 
    uniqueReal() const;

    // Returns the IndexType
    IndexType 
    type() const;

    // Returns the name of this Index
    std::string 
    name() const;

    // Returns the name of this Index with primes removed
    const std::string&
    rawname() const;

    // Returns true if Index default initialized
    bool 
    isNull() const;

    // Returns the Arrow direction of this Index
    Arrow 
    dir() const { return Out; }

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
    operator<(const Index& other) const;

    // Creates an IndexVal from this Index with value i
    IndexVal 
    operator()(int i) const;

    //
    // Other methods
    //

    // Conjugate this Index.
    // Currently has no effect; exists for forward compatibility
    // with Arrows and interface compatibility with class IQIndex.
    void 
    conj() { } //for forward compatibility with arrows

    // Write Index to binary output stream.
    void 
    write(std::ostream& s) const;

    // Read Index from binary input stream.
    Index& 
    read(std::istream& s);

    //
    // Static Index instances
    //

    //Static default-constructed placeholder Index
    static const 
    Index& Null();

    private:

    /////////////
    IndexDatPtr p;

    int primelevel_; 
    /////////////

    explicit
    Index(const IndexDatPtr& p, int plev);

    }; //class Index


//
// IndexVal
//
// Class pairing an Index of dimension m
// with a specific value i where 1 <= i <= m
//
class IndexVal : public Index
    {
    public:

    int i;

    IndexVal();

    IndexVal(const Index& index, int i_);

    bool
    operator==(const IndexVal& other) const
        {
        return (Index::operator==(other) && i == other.i);
        }

    bool
    operator!=(const IndexVal& other) const
        {
        return !operator==(other);
        }

    static const IndexVal& 
    Null();
    };


//Return a copy of I, increasing primelevel.
template<class T>
T
primed(T I, int inc = 1) { I.prime(inc); return I; }

//Return a copy of I, increasing primelevel if I.type() == type
template<class T>
T 
primed(T I, IndexType type, int inc = 1) { I.prime(type,inc); return I; }

//Return a copy of I with primelevel set to zero.
template<class T>
T
deprimed(T I, IndexType type = All) { I.noprime(type); return I; }

//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
template<class T>
T
mapPrime(T I, int plevold, int plevnew, IndexType type = All)
    { I.mapprime(plevold,plevnew,type); return I; }

template <class T>
T
conj(T res) { res.conj(); return res; }

//Returns a string version of this Index's bond dimension.
std::string
showm(const Index& I);

std::string 
nameint(const std::string& f, int n);

std::ostream& 
operator<<(std::ostream & s, const Index &t);

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv);

std::ostream& 
operator<<(std::ostream& s, const IndexType& it);

#undef Cout
#undef Format
#undef Endl

#endif
