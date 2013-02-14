//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEX_H
#define __ITENSOR_INDEX_H
#include <string>
#include "global.h"
#include "boost/shared_ptr.hpp"
#include "boost/uuid/uuid.hpp"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

enum IndexType { Link, Site, ReIm, All };

//Forward declarations
class IndexDat;
struct IndexVal;


//
// Index
//
// Represents a tensor index of fixed bond dimension m.
// The == operator can be used to determine if two Index's match
// (are copies of the same Index instance).
// An Index can be made temporarily distinct from other copies 
// by increasing its primeLevel.
//
class Index
    {
    public:
    //
    // Constructors
    //

    // For a default constructed Index, isNull() returns true.
    Index();

    // Standard Index constructor.
    // Name of Index is used for printing purposes
    Index(const std::string& name, 
          int m = 1, 
          IndexType it = Link, 
          int primelev = 0);

    // Copy constructor which increments the primelevel of the copy.
    Index(IndexType pt,
          const Index& other, 
          int primeinc = 1);

    //
    // Accessor Methods
    //

    // Returns the bond dimension of an Index.
    int 
    m() const;

    // Returns the IndexType of this Index.
    IndexType 
    type() const;

    // Returns the name of this Index.
    std::string 
    name() const;

    // Returns the name of this Index with primes removed.
    const std::string&
    rawname() const;

    // Returns a unique Real number identifying this Index.
    // Useful for rapidly checking that two Index instances match.
    Real 
    uniqueReal() const;

    // Returns true if Index was default initialized.
    bool 
    isNull() const;
    // Returns true if Index was NOT default initialized.
    bool 
    isNotNull() const;

    // Returns the prime level of this Index
    int 
    primeLevel() const;
    // Set the prime level to a specified value.
    void 
    primeLevel(int plev);

    // Returns the Arrow direction of this Index.
    Arrow 
    dir() const { return Out; }

    //
    // Operators
    //

    // Equality (==) operator.
    // True if other Index is a copy of this Index with the same primelevel.
    bool 
    operator==(const Index& other) const;

    bool 
    operator!=(const Index& other) const { return !operator==(other); }

    // Check if other Index is a copy of this, ignoring primeLevel.
    bool 
    noprimeEquals(const Index& other) const;

    // Useful for sorting Index objects.
    bool 
    operator<(const Index& other) const;

    // Creates an IndexVal from this Index with index i.
    IndexVal 
    operator()(int i) const;

    friend std::ostream& 
    operator<<(std::ostream & s, const Index &t);

    //
    // Prime methods
    //

    // Increase primelevel by 1 (or by optional amount inc).
    void 
    prime(int inc = 1);

    // Increase primelevel by 1 (or optional amount inc)
    // if type matches this Index or type==All
    void 
    prime(IndexType type, int inc = 1);

    // Set primelevel to zero (optionally only if type matches)
    void 
    noprime(IndexType type = All) { prime(type,-primelevel_); }

    // Switch primelevel from plevold to plevnew. 
    // Has no effect if plevold doesn't match current primelevel.
    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //
    // Other methods
    //

    // Write Index to binary output stream.
    void 
    write(std::ostream& s) const;

    // Read Index from binary input stream.
    void 
    read(std::istream& s);

    // Print Index to stdout.
    void 
    print(const std::string& name = "") const
        { std::cout << "\n" << name << " =\n" << *this << std::endl; }

    // Conjugate this Index.
    // Currently has no effect; for forward compatibility
    // with Arrows and interface compatibility with class IQIndex.
    void 
    conj() { } //for forward compatibility with arrows


    static const Index& Null();

    // Static Index indexing real and imaginary parts of a complex ITensor.
    static const Index& IndReIm();

    // IndReIm with primeLevel 1
    static const Index& IndReImP();

    // IndReIm with primeLevel 2
    static const Index& IndReImPP();

    enum Imaker { makeReIm, makeReImP, makeReImPP, makeNull };

    private:

    friend class IQIndex;

    explicit
    Index(Imaker im);

    /////////////
    //
    // Data Members

    boost::shared_ptr<IndexDat> p;

    int primelevel_; 

    //
    /////////////

    }; //class Index


//
// IndexVal
//
// Struct pairing an Index (of dimension m)
// with a specific value i where 1 <= i <= m.
//
struct IndexVal
    {
    Index ind; 
    int i;

    IndexVal();

    IndexVal(const Index& index, int i_);

    bool 
    operator==(const IndexVal& other) const; 

    void
    prime(int inc = 1) { ind.prime(inc); }

    void
    prime(IndexType type, int inc = 1) { ind.prime(type,inc); }

    static const IndexVal& 
    Null();

    private:

    explicit
    IndexVal(Index::Imaker im);

    };


// Make a copy of this Index, increasing primelevel.
Index inline
primed(Index I, int inc = 1) { I.prime(inc); return I; }

Index inline
primed(Index I, IndexType type, int inc = 1) { I.prime(type,inc); return I; }

// Return a copy of this Index with primelevel set to zero.
Index inline
deprimed(Index I) { I.noprime(); return I; }

// Returns a string version of this Index's bond dimension.
std::string
showm(const Index& I);

template <class T> T 
conj(T res) 
    { 
    res.conj(); 
    return res; 
    }

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv);

std::ostream& 
operator<<(std::ostream& s, const boost::uuids::uuid& id);

std::ostream& 
operator<<(std::ostream& s, const IndexType& it);

std::string 
nameint(const std::string& f, int n);

#undef Cout
#undef Format
#undef Endl

#endif
