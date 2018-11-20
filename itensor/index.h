//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEX_H
#define __ITENSOR_INDEX_H
#include "itensor/global.h"
#include "itensor/smallstring.h"
#include "itensor/tagset.h"
#include "itensor/arrow.h"
#include <thread>

namespace itensor {

//Forward declarations
class IndexVal;

namespace detail {
    struct RandomID
        {
        using rng_type = std::mt19937_64;
        using result_type = typename rng_type::result_type;
        std::hash<std::thread::id> hasher;
        RandomID()
            : rng(std::clock() + hasher(std::this_thread::get_id()))
            { }

        result_type
        operator()() { return rng(); }

        private:
        rng_type rng;
        };

    struct SequentialID
        {
        using result_type = std::uint_fast32_t;

        SequentialID() : id(1) { }

        result_type
        operator()()
            { 
            auto res = id;
            id += 1;
            return res; 
            }

        private:
        result_type id = 0;
        };
    } //namespace detail

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
    using IDGenerator = detail::RandomID;
    using id_type = IDGenerator::result_type;
    using indexval_type = IndexVal;
    using prime_type = int;
    using extent_type = int;
    private:
    id_type id_;
    extent_type m_;
    prime_type primelevel_; 
    TagSet tags_;
    public:

    Index();

    explicit
    Index(long m, 
          const TagSet& t = "",
          int primelev = 0);

    explicit
    Index(long m, 
          int primelev);

    // Returns the bond dimension
    long 
    m() const { return m_; }

    // Returns the prime level
    int 
    primeLevel() const { return primelevel_; }
    // Sets the prime level to a specified value.
    Index& 
    primeLevel(int plev);

    // Returns the TagSet
    TagSet
    tags() const { return tags_; }

    id_type
    id() const { return id_; }

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

    Index& 
    setPrime(int p) { primelevel_ = p; return *this; };

    Index& 
    noPrime() { primelevel_ = 0; return *this; };

    // Increase primelevel by 1 (or by optional amount inc)
    Index& 
    prime(int inc = 1);

    //TODO: clean up
    // Increase primelevel by 1 (or optional amount inc)
    // if type matches this Index or type==All
    //Index& 
    //prime(IndexType type, int inc = 1);

    // Set primelevel to zero (optionally only if type matches)
    //Index& 
    //noprime(IndexType type = All) { prime(type,-primelevel_); return *this; }

    // Switch primelevel from plevold to plevnew
    // Has no effect if plevold doesn't match current primelevel
    //Index& 
    //mapprime(int plevold, int plevnew, IndexType type = All);

    // Check if other Index is a copy of this, ignoring primeLevel
    bool 
    noprimeEquals(Index const& other) const;

    // Add tags
    Index&
    addTags(const TagSet& t) { tags_.addTags(t); return *this; }

    // Remove tags
    Index&
    removeTags(const TagSet& t) { tags_.removeTags(t); return *this; }

    // Set tags
    Index&
    setTags(const TagSet& t) { tags_ = t; return *this; }

    // Return a copy of this Index with new tags
    Index
    operator()(const TagSet& t) const { auto I = *this; I.setTags(t); return I; }

    //Return an IndexVal with specified value
    IndexVal
    operator()(long val) const;

    //Return copy of this Index with primelevel plev
    Index
    operator[](int plev) const;

    // Conjugate this Index.
    // Currently has no effect; exists for forward compatibility
    // with Arrows and interface compatibility with class IQIndex.
    void 
    dag() { } //for forward compatibility with arrows

    //define size()==m() in order to do 
    //for(auto n : range(I)) { ... } for some Index I
    long
    size() const { return m(); }

    // Write Index to binary output stream.
    void 
    write(std::ostream& s) const;

    // Read Index from binary input stream.
    Index& 
    read(std::istream& s);

    private:

    Index::id_type 
    generateID();

    }; //class Index

// i1 compares equal to i2 if i2 is a copy of i1 with same primelevel
bool 
operator==(Index const& i1, Index const& i2);
bool 
operator!=(Index const& i1, Index const& i2);

// Useful for sorting Index objects
bool 
operator<(Index const& i1, Index const& i2);
bool 
operator>(Index const& i1, Index const& i2);


//
// IndexVal
//
// Class pairing an Index of dimension m
// with a specific value i where 1 <= i <= m
//
class IndexVal
    {
    public:
    using index_type = Index;

    Index index;
    long val;

    IndexVal();

    IndexVal(const Index& index, long val_);

    long
    m() const { return index.m(); }

    explicit operator bool() const { return bool(index); }

    IndexVal& 
    prime(int inc = 1) { index.prime(inc); return *this; }

    //TODO: clean up
    //IndexVal& 
    //prime(IndexType type, int inc = 1) { index.prime(type,inc); return *this; }

    IndexVal& 
    noPrime() { index.noPrime(); return *this; }

    //TODO: clean up
    //IndexVal& 
    //mapprime(int plevold, int plevnew, IndexType type = All) 
    //    { index.mapprime(plevold,plevnew,type); return *this; }

    void
    dag() { }

    };

bool
operator==(IndexVal const& iv1, IndexVal const& iv2);
bool
operator!=(IndexVal const& iv1, IndexVal const& iv2);

bool
operator==(IndexVal const& iv, Index const& I);
bool
operator==(Index const& I, IndexVal const& iv);

TagSet inline
tags(const Index& I) { return I.tags(); }

Index inline
addTags(Index I, const TagSet& t) { I.addTags(t); return I; }

Index inline
removeTags(Index I, const TagSet& t) { I.removeTags(t); return I; }

Index inline
setTags(Index I, const TagSet& t) { I.setTags(t); return I; }

//
// Check if Index I contains the tags tsmatch.
// If tsmatch==TagSet(All), return true
//
bool inline
hasTags(Index I, const TagSet& tsmatch) { return tsmatch==TagSet(All) || hasTags(tags(I),tsmatch); }

//
// Return true if Index I contains tags tsmatch and has prime level plmatch
// If tsmatch==TagSet(All), Index I can have any tags
// If plmatch < 0, Index I can have any prime level
//
bool inline
matchTagsPrime(Index I, TagSet const& tsmatch, int plmatch) { return hasTags(I,tsmatch) && (plmatch<0 || plmatch==I.primeLevel()); }

//
// Return true if Index I has tags tsmatch and has prime level plmatch
// If tsmatch==TagSet(All), Index I can have any tags
// If plmatch < 0, Index I can have any prime level
//
bool inline
matchTagsPrimeExact(Index I, TagSet const& tsmatch, int plmatch) { return tags(I)==tsmatch && (plmatch<0 || plmatch==I.primeLevel()); }

Index inline
dag(Index res) { res.dag(); return res; }

IndexVal inline
dag(IndexVal res) { res.dag(); return res; }

template<typename... VarArgs>
Index
prime(Index I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
Index
noPrime(Index I, VarArgs&&... vargs) { I.noPrime(std::forward<VarArgs>(vargs)...); return I; }

//TODO: clean up
//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
//Index inline
//mapprime(Index I, int plevold, int plevnew, IndexType type = All)
//    { I.mapprime(plevold,plevnew,type); return I; }

template<typename... VarArgs>
IndexVal
prime(IndexVal I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
IndexVal
noPrime(IndexVal I, VarArgs&&... vargs) { I.noPrime(std::forward<VarArgs>(vargs)...); return I; }

//TODO: clean up
//Return a copy of I with prime level changed to plevnew if
//old prime level was plevold. Otherwise has no effect.
//IndexVal inline
//mapprime(IndexVal I, int plevold, int plevnew, IndexType type = All)
//    { I.mapprime(plevold,plevnew,type); return I; }

//Make a new index with same properties as I,
//but a different id number (will not compare equal)
//and primelevel zero (or specified value)
Index
sim(Index const& I, int plev = 0);

//Returns a string version of this Index's bond dimension.
std::string
showm(Index const& I);

//TODO: clean up
//Depecreate, nameint is a strange name when Indices don't
//have names anymore, easy enough to write format("%s%d",f,d)
//std::string 
//nameint(std::string const& f, int n);

std::ostream& 
operator<<(std::ostream & s, Index const& t);

std::ostream& 
operator<<(std::ostream& s, IndexVal const& iv);

void
add(Args& args, 
    Args::Name const& name, 
    TagSet const& ts);

TagSet
getTagSet(Args const& args, 
          Args::Name const& name);

TagSet
getTagSet(Args const& args, 
          Args::Name const& name, 
          TagSet const& default_val);

} //namespace itensor

#endif
