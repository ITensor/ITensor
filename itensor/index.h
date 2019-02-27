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
#include "itensor/qn.h"
#include <thread>

namespace itensor {

//Forward declarations
class IndexVal;

using QNInt = std::pair<QN,long>;

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

class IQIndexDat;

class Index
    {
    public:
    using IDGenerator = detail::RandomID;
    using id_type = IDGenerator::result_type;
    using indexval_type = IndexVal;
    using prime_type = int;
    using extent_type = int;

    using qnstorage = std::vector<QNInt>;
    using qn_ptr = std::shared_ptr<IQIndexDat>;

    private:
    id_type id_;
    extent_type dim_;
    prime_type primelevel_; 
    Arrow dir_ = Out;
    qn_ptr pd;
    TagSet tags_;
    public:

    Index();

    explicit
    Index(long m, 
          TagSet const& ts = TagSet());

    explicit
    Index(std::string s,
          long m)
        {
        Error("Index(string,int) constructor deprecated, use Index(int,...) instead");
        }

    template<typename... QN_Sizes>
    Index(QN const& q1, long size1,
          QN_Sizes const&... qnsizes);

    Index(qnstorage && qns, 
          TagSet const& ts = TagSet());

    Index(qnstorage && qns, 
          Arrow dir,
          TagSet const& ts);

    // Returns the dimension of this Index
    long 
    dim() const { return dim_; }

    // Returns the prime level
    int 
    primeLevel() const { return primelevel_; }

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

    // Sets the prime level to a specified value.
    Index& 
    setPrime(int p);

    // Sets the prime level to 0.
    Index& 
    noPrime();

    // Increase primelevel by 1 (or by optional amount inc)
    Index& 
    prime(int inc = 1);

    // Add tags
    Index&
    addTags(const TagSet& t) { tags_.addTags(t); return *this; }

    // Remove tags
    Index&
    removeTags(const TagSet& t) { tags_.removeTags(t); return *this; }

    // Set tags
    Index&
    setTags(const TagSet& t) { tags_ = t; return *this; }

    // Set tags
    Index&
    replaceTags(const TagSet& tsold, const TagSet& tsnew) { tags_.removeTags(tsold); tags_.addTags(tsnew); return *this; }

    // Return a copy of this Index with new tags
    Index
    operator()(const TagSet& t) const { auto I = *this; I.setTags(t); return I; }

    //Return an IndexVal with specified value
    IndexVal
    operator()(long val) const;
    IndexVal
    operator=(long val) const;

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

    //
    // QN related functions
    // 
    
    //number of quantum number blocks
    long 
    nblock() const;

    //1-indexed
    long
    blocksize(long i) const;
      
    //1-indexed
    QN const& 
    qn(long i) const;

    Arrow 
    dir() const { return dir_; }
    void
    dir(Arrow ndir) { dir_ = ndir; }

    Index& 
    dag() { dir_ = -dir_; return *this; }

    qn_ptr const&
    store() const { return pd; }

    private:

    void
    makeStorage(qnstorage && qi);

    Index::id_type 
    generateID();

    public:

    //
    // Advanced / developer methods.
    // Not intended for normal usage.
    //

    // Constructor taking a QN pointer
    Index(qn_ptr const& p,
          TagSet const& tags = TagSet());

    Index(qn_ptr const& p,
          Arrow dir, 
          TagSet const& tags);

    //0-indexed
    long
    blocksize0(long i) const;

    void
    removeQNs() { pd.reset(); }

    public:

    //Deprecated: prefer to use I.dim() or dim(I)
    long 
    m() const { return dim_; }

    }; //class Index

// i1 compares equal to i2 if i2 is a copy of i1 with same primelevel
bool 
operator==(Index const& i1, Index const& i2);
bool 
operator!=(Index const& i1, Index const& i2);

bool
equalsIgnorePrime(Index const& i1, Index const& i2);

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
    dim() const { return index.dim(); }

    explicit operator bool() const { return bool(index); }

    IndexVal& 
    prime(int inc = 1) { index.prime(inc); return *this; }

    IndexVal&
    setPrime(int n) { index.setPrime(n); return *this; }

    IndexVal& 
    noPrime() { index.noPrime(); return *this; }

    IndexVal& 
    dag();

    QN const&
    qn() const;


    //Deprecated: prefer to use .dim()
    long
    m() const { return index.m(); }

    };

bool
operator==(IndexVal const& iv1, IndexVal const& iv2);
bool
operator!=(IndexVal const& iv1, IndexVal const& iv2);

bool
operator==(IndexVal const& iv, Index const& I);
bool
operator==(Index const& I, IndexVal const& iv);

long inline
dim(Index const& I) { return I.dim(); }

long inline
dim(IndexVal const& I) { return I.dim(); }

int inline
primeLevel(Index const& I) { return I.primeLevel(); }

TagSet inline
tags(const Index& I) { return I.tags(); }

Index inline
addTags(Index I, const TagSet& t) { I.addTags(t); return I; }

Index inline
removeTags(Index I, const TagSet& t) { I.removeTags(t); return I; }

Index inline
setTags(Index I, const TagSet& t) { I.setTags(t); return I; }

Index inline
replaceTags(Index I, const TagSet& tsold, const TagSet& tsnew) { I.replaceTags(tsold,tsnew); return I; }

Index
tags(Index I, std::string st);

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
matchTagsPrime(Index I, TagSet const& tsmatch, int plmatch) { return hasTags(I,tsmatch) && (plmatch<0 || plmatch==primeLevel(I)); }

//
// Return true if Index I has tags tsmatch and has prime level plmatch
// If tsmatch==TagSet(All), Index I can have any tags
// If plmatch < 0, Index I can have any prime level
//
bool inline
matchTagsPrimeExact(Index I, TagSet const& tsmatch, int plmatch) { return tags(I)==tsmatch && (plmatch<0 || plmatch==primeLevel(I)); }

bool inline
hasQNs(Index const& I) { return I.nblock()!=0; }

Index inline
removeQNs(Index I) { if(hasQNs(I)) I.removeQNs(); return I; }

bool inline
hasQNs(IndexVal const& iv) { return hasQNs(iv.index); }
  
Index inline
dag(Index res) { res.dag(); return res; }

IndexVal inline
dag(IndexVal res) { res.dag(); return res; }

template<typename... VarArgs>
Index
prime(Index I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
Index
setPrime(Index I, VarArgs&&... vargs) { I.setPrime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
Index
noPrime(Index I, VarArgs&&... vargs) { I.noPrime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
IndexVal
prime(IndexVal I, VarArgs&&... vargs) { I.prime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
IndexVal
setPrime(IndexVal I, VarArgs&&... vargs) { I.setPrime(std::forward<VarArgs>(vargs)...); return I; }

template<typename... VarArgs>
IndexVal
noPrime(IndexVal I, VarArgs&&... vargs) { I.noPrime(std::forward<VarArgs>(vargs)...); return I; }


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

long
QNblock(Index const& I, 
        QN const& Q);

long
QNblockSize(Index const& I, 
            QN const& Q);

void
write(std::ostream & s, QNInt const& q);

void
read(std::istream & s, QNInt & q);

} //namespace itensor

#include "itensor/index_impl.h"

#endif
