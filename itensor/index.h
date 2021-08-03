//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_INDEX_H
#define __ITENSOR_INDEX_H
#include "itensor/global.h"
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
    using extent_type = int;

    using qnstorage = std::vector<QNInt>;
    using qn_ptr = std::shared_ptr<IQIndexDat>;

    private:
    id_type id_;
    extent_type dim_;
    Arrow dir_ = Out;
    qn_ptr pd;
    TagSet tags_;
    public:

    Index();

    explicit
    Index(long dim, 
          TagSet const& ts = TagSet("0"));

    // Deprecated
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
          TagSet const& ts = TagSet("0"));

    Index(qnstorage && qns, 
          Arrow dir,
          TagSet const& ts = TagSet("0"));

    // Returns the dimension of this Index
    long 
    dim() const { return dim_; }

    // Returns the prime level
    int 
    primeLevel() const { return tags_.primeLevel(); }

    // Returns the TagSet
    TagSet
    tags() const { return tags_; }

    id_type
    id() const { return id_; }

    // Evaluates to false if Index is default constructed.
    explicit operator bool() const;

    // (Explicitly) convertible to integer types
    explicit operator int() const { return dim(); }
    explicit operator long() const { return dim(); }
    explicit operator size_t() const { return dim(); }

    // Add tags
    Index&
    addTags(const TagSet& t) { tags_.addTags(t); return *this; }

    // Remove tags
    Index&
    removeTags(const TagSet& t) { tags_.removeTags(t); return *this; }

    // Set tags
    Index&
    setTags(const TagSet& t) { tags_.setTags(t); return *this; }

    // Remove all tags
    Index&
    noTags() { tags_.noTags(); return *this; }

    // Set tags
    Index&
    replaceTags(const TagSet& tsold, const TagSet& tsnew) { tags_.replaceTags(tsold,tsnew); return *this; }

    // Sets the prime level to a specified value.
    Index& 
    setPrime(int p);

    // Sets the prime level to 0.
    Index& 
    noPrime();

    // Increase primelevel by 1 (or by optional amount inc)
    Index& 
    prime(int inc = 1);

    //Return an IndexVal with specified value
    IndexVal
    operator()(long val) const;
    IndexVal
    operator=(long val) const;

    //define size()==dim() in order to do 
    //for(auto n : range(I)) { ... } for some Index I
    long
    size() const { return dim(); }

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
    setDir(Arrow ndir) { dir_ = ndir; }

    Index& 
    dag() { dir_ = -dir_; return *this; }

    qn_ptr const&
    store() const { return pd; }

    Index&
    sim(Index const& I)
      {
      *this = I;
      id_ = generateID();
      return *this;
      }

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

    Index(id_type id,
          long dim, 
          Arrow dir, 
          TagSet const& ts);

    Index(id_type id,
          long dim, 
          Arrow dir, 
          TagSet const& ts,
          qnstorage&& qns);

    //0-indexed
    long
    blocksize0(long i) const;

    void
    removeQNs() { pd.reset(); }

    public:

    //Deprecated: prefer to use I.dim() or dim(I)
    //long 
    //m() const
    //    {
    //    Global::warnDeprecated(".m() is deprecated in favor of dim(Index)");
    //    return this->dim();
    //    }

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

    Index index;
    long val;

    IndexVal();

    IndexVal(const Index& index, long val_);

    long
    dim() const { return index.dim(); }

    explicit operator bool() const { return bool(index); }

    // Add tags
    IndexVal&
    addTags(const TagSet& t) { index.addTags(t); return *this; }

    // Remove tags
    IndexVal&
    removeTags(const TagSet& t) { index.removeTags(t); return *this; }

    // Set tags
    IndexVal&
    setTags(const TagSet& t) { index.setTags(t); return *this; }

    // Remove all tags
    IndexVal&
    noTags() { index.noTags(); return *this; }

    // Set tags
    IndexVal&
    replaceTags(const TagSet& tsold, const TagSet& tsnew) { index.replaceTags(tsold,tsnew); return *this; }

    // Sets the prime level to a specified value.
    IndexVal&
    setPrime(int p) { index.setPrime(p); return *this; }

    // Sets the prime level to 0.
    IndexVal&
    noPrime() { index.noPrime(); return *this; }

    // Increase primelevel by 1 (or by optional amount inc)
    IndexVal&
    prime(int inc = 1) { index.prime(inc); return *this; }

    IndexVal& 
    dag();

    QN const&
    qn() const;


    //Deprecated: prefer to use .dim()
    long
    m() const 
      {
      Global::warnDeprecated(".m() is deprecated in favor of dim(IndexVal)");
      return this->dim();
      }

    };

bool
operator==(IndexVal const& iv1, IndexVal const& iv2);
bool
operator!=(IndexVal const& iv1, IndexVal const& iv2);

bool
operator==(IndexVal const& iv, Index const& I);
bool
operator==(Index const& I, IndexVal const& iv);

Index::id_type
id(Index const& I);

long
dim(Index const& I);

long
dim(IndexVal const& I);

int
primeLevel(Index const& I);

TagSet
tags(Index const& I);

//
// Check if Index I contains the tags tsmatch.
//
bool
hasTags(Index I, TagSet const& tsmatch);

bool
hasQNs(Index const& I);

Index
removeQNs(Index I);

bool
hasQNs(IndexVal const& iv);
  
Index
dag(Index res);

IndexVal
dag(IndexVal res);

QN const&
qn(Index const& i, long b);

QN const&
qn(IndexVal const& iv);

Arrow
dir(Index const& res);

Arrow
dir(IndexVal const& res);

Index const&
index(IndexVal const& res);

long
val(IndexVal const& res);

long
nblock(Index const& i);

long
blocksize(Index const& i, long b);

//
// Tag functions
//

Index
addTags(Index I, TagSet const& t);

Index
removeTags(Index I, TagSet const& t);

Index
setTags(Index I, TagSet const& t);

Index
noTags(Index I);

Index
replaceTags(Index I, TagSet const& tsold, TagSet const& tsnew);

Index
prime(Index I, int plinc = 1);

Index
setPrime(Index I, int plev);

Index
noPrime(Index I);

IndexVal
addTags(IndexVal I, TagSet const& t);

IndexVal
removeTags(IndexVal I, TagSet const& t);

IndexVal
setTags(IndexVal I, TagSet const& t);

IndexVal
noTags(IndexVal I);

IndexVal
replaceTags(IndexVal I, TagSet const& tsold, TagSet const& tsnew);

IndexVal
prime(IndexVal I, int plinc = 1); 

IndexVal
setPrime(IndexVal I, int plev); 

IndexVal
noPrime(IndexVal I); 

// Get the direct sum of two indices
Index
directSum(Index const& i,
          Index const& j,
          Args const& args = Args::global());

//Make a new index with same properties as I,
//but a different id number (will not compare equal)
Index
sim(Index const& I);

//Returns a string version of this Index's bond dimension.
std::string
showDim(Index const& I);

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

bool
isFermionic(Index const& I);


#ifdef ITENSOR_USE_HDF5
void
h5_write(h5::group parent, std::string const& name, Index const& I);
void
h5_read(h5::group parent, std::string const& name, Index & I);
#endif

} //namespace itensor

#include "itensor/index_impl.h"

#endif
