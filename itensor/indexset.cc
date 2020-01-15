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
#include "itensor/indexset.h"

namespace itensor {

void
checkIndexSet(IndexSet const& is)
    {
    //Check if any duplicate indices
    for(size_t j = 0; j < is.size(); ++j)
    for(size_t k = j+1; k < is.size(); ++k)
        if( is[j] == is[k] )
            {
            println("index set = \n",is);
            throw ITError("Duplicate indices in index set");
            }
    }

void
checkIndexPositions(std::vector<int> const& is)
    {
    //Check if there any duplicate index positions
    //(it is ok if Indices are not found, they are just ignored)
    for(size_t j = 0; j < is.size(); ++j)
    for(size_t k = j+1; k < is.size(); ++k)
        if( (is[j] != -1) && (is[j] == is[k]) )
            throw ITError("An Index was found more than once in the IndexSet");
    }

IndexSet& IndexSet::
dag() { for(auto& J : *this) J.dag(); return *this; }

IndexSet
dag(IndexSet is)
  {
  is.dag();
  return is;
  }

long
order(IndexSet const& is)
    {
    return is.order();
    }

long
length(IndexSet const& is)
    {
    return order(is);
    }

IndexSet::iterator IndexSet::
begin() { return iterator{*this}; }

IndexSet::iterator IndexSet::
end() { return iterator::makeEnd(*this); }

IndexSet::const_iterator IndexSet::
begin() const { return const_iterator{*this}; }

IndexSet::const_iterator IndexSet::
end() const { return const_iterator::makeEnd(*this); }

IndexSet::const_iterator IndexSet::
cbegin() const { return begin(); }

IndexSet::const_iterator IndexSet::
cend() const { return end(); }

IndexSet& IndexSet::
setTags(TagSet const& tsnew)
    {
    for(auto& J : *this) J.setTags(tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
setTags(TagSet const& tsnew,
        IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).setTags(tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
setTags(TagSet const& tsnew,
        TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.setTags(tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
noTags()
    {
    for(auto& J : *this) J.noTags();
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
noTags(IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).noTags();
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
noTags(TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.noTags();
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
addTags(TagSet const& tsadd)
    {
    for(auto& J : *this) J.addTags(tsadd);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
addTags(TagSet const& tsadd,
        IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).addTags(tsadd);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
addTags(TagSet const& tsadd,
        TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.addTags(tsadd);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
removeTags(TagSet const& tsremove)
    {
    for(auto& J : *this) J.removeTags(tsremove);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
removeTags(TagSet const& tsremove,
           IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).removeTags(tsremove);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
removeTags(TagSet const& tsremove,
           TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.removeTags(tsremove);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew)
    {
    for(auto& J : *this) J.replaceTags(tsold,tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew,
            IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).replaceTags(tsold,tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew,
            TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.replaceTags(tsold,tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
setPrime(int plnew)
    {
    for(auto& J : *this) J.setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
setPrime(int plnew, IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
setPrime(int plnew, TagSet const& tsmatch)
    {
    for(auto& J : *this) if( hasTags(J,tsmatch) ) J.setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
mapPrime(int plold, int plnew)
    {
    for(auto& J : *this) if( J.primeLevel() == plold ) J.setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
mapPrime(int plold, int plnew, IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs)
      {
      auto& J = parent::index(i);
      if( J.primeLevel() == plold ) J.setPrime(plnew);
      }
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
mapPrime(int plold, int plnew, TagSet const& tsmatch)
    {
    for(auto& J : *this)
      if( (J.primeLevel() == plold) && hasTags(J,tsmatch) ) J.setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
prime(int plinc)
    {
    for(auto& J : *this) J.prime(plinc);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
prime(int plinc, IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).prime(plinc);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet& IndexSet::
prime(int plinc, TagSet const& tsmatch)
    {
    for(auto& J : *this) if( hasTags(J,tsmatch) ) J.prime(plinc);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    return *this;
    }

IndexSet
sim(IndexSet is)
    {
    for(auto& J : is)
        J = sim(J);
    return is;
    }

IndexSet& IndexSet::
removeQNs()
    {
    for(auto& J : *this) J.removeQNs();
    return *this;
    }

//
// Methods for Manipulating IndexSet
//

void
write(std::ostream& s, IndexSet const& is)
    {
    using parent = typename IndexSet::parent;
    parent const& pr = is;
    itensor::write(s,pr);
    }

void
read(std::istream& s, IndexSet & is)
    {
    using parent = typename IndexSet::parent;
    parent & pr = is;
    itensor::read(s,pr);
    }

Arrow
dir(IndexSet const& is, Index const& I)
    {
    for(auto const& J : is)
        {
        if(J == I) return J.dir();
        }
    throw ITError("dir: Index not found");
    return In;
    }


//
// Given IndexSet iset and Index I,
// return int j such that iset[j] == I.
// If not found, returns -1
//
int
indexPosition(IndexSet const& is,
              Index const& imatch)
    {
    for(auto j : range(order(is)))
        if( is[j] == imatch ) return j;
    return -1;
    }


std::vector<int>
indexPositions(IndexSet const& is,
               IndexSet const& ismatch)
    {
    auto ilocs = std::vector<int>();
    for(auto const& J : ismatch)
        {
        auto loc = indexPosition(is,J);
        if( loc != -1 ) ilocs.push_back(loc);
        }
#ifdef DEBUG
    checkIndexPositions(ilocs);
#endif
    return ilocs;
    }

bool
hasIndex(IndexSet const& is,
         Index const& imatch)
  {
  for( auto& I : is )
    if( I == imatch ) return true;
  return false;
  }

bool
hasInds(IndexSet const& is,   
        IndexSet const& ismatch)
  {
  for( auto& I : ismatch )
    if( !hasIndex(is,I) ) return false;
  return true;
  }

bool
hasSameInds(IndexSet const& is1,
            IndexSet const& is2)
  {
  return (order(is1)==order(is2)) && hasInds(is1,is2);
  }

bool
equals(IndexSet const& is1,
       IndexSet const& is2)
  {
  auto N = order(is1);
  if( order(is2)!= N ) return false;
  for( auto n : range1(N) )
    if( is1(n)!=is2(n) ) return false;
  return true;
  }

long
minDim(IndexSet const& iset)
    {
    if(iset.empty()) return 1l;
    auto mm = dim(iset[0]);
    for(long j = 1; j < iset.order(); ++j)
        mm = std::min(mm,dim(iset[j]));

    return mm;
    }

long
maxDim(IndexSet const& iset)
    {
    if(iset.empty()) return 1l;

    auto mm = dim(iset[0]);
    for(long j = 1; j < iset.order(); ++j)
        mm = std::max(mm,dim(iset[j]));

    return mm;
    }

IndexSet
sim(IndexSet is,
    IndexSet const& ismatch)
    {
    for(auto& J : is)
        if( hasIndex(ismatch,J) )
            J = sim(J);
    return is;
    }

IndexSet
sim(IndexSet is,
    TagSet const& tsmatch)
    {
    for(auto& J : is)
        if( hasTags(J,tsmatch) )
            J = sim(J);
    return is;
    }

void
contractIS(IndexSet const& Lis,
           IndexSet const& Ris,
           IndexSet & Nis,
           bool sortResult)
    {
    Labels Lind,
           Rind,
           Nind;
    computeLabels(Lis,Lis.order(),Ris,Ris.order(),Lind,Rind);
    contractIS(Lis,Lind,Ris,Rind,Nis,Nind,sortResult);
    }

std::ostream&
operator<<(std::ostream& s, IndexSet const& is)
    {
    for(auto i : range1(is.order()))
        {
        if(hasQNs(is)) s << is(i) << std::endl;
        else s << is(i) << " ";
        }
    return s;
    }

bool
hasQNs(IndexSet const& is)
    {
    for(auto& I : is) if(hasQNs(I)) return true;
    return false;
    }

void
checkQNConsistent(IndexSet const& is)
    {
#ifdef DEBUG
    if(is.order() > 0 && hasQNs(is.front()))
        {
        for(long n = 1; n < is.order(); n += 1)
            {
            if(not hasQNs(is[n]))
                {
                println("Non-QN index = ",is[n]);
                Error("IndexSet: cannot mix QN and non-QN Indices");
                }
            }
        }
#endif
    }

IndexSet
findInds(IndexSet const& is,
         TagSet const& tsmatch)
    {
    auto inds = std::vector<Index>();
    for( auto& I : is )
        if( hasTags(I,tsmatch) )
            inds.push_back(std::move(I));
    return IndexSet(inds);
    }

Index
findIndex(IndexSet const& is)
    {
    if( order(is) > 1 ) Error("In findIndex: more than one Index found, consider using findInds instead");
    else if( order(is) == 1 ) return is(1);
    return Index();
    }

Index
findIndex(IndexSet const& is,
          TagSet const& tsmatch)
    {
    return findIndex(findInds(is,tsmatch));
    }

IndexSet
findIndsExcept(IndexSet const& is,
               TagSet const& tsmatch)
    {
    auto inds = std::vector<Index>();
    for( auto& I : is )
        if( !hasTags(I,tsmatch) )
            inds.push_back(std::move(I));
    return IndexSet(inds);
    }

IndexSet
commonInds(IndexSet const& is1,
           IndexSet const& is2)
    {
    auto inds = std::vector<Index>();
    for( auto& I : is1 )
        if( hasIndex(is2,I) )
            inds.push_back(std::move(I));
    return IndexSet(inds);
    }

IndexSet
unionInds(IndexSet const& is1,
          IndexSet const& is2)
    {
    auto inds = std::vector<Index>();
    for( auto& I : is1 )
        inds.push_back(std::move(I));
    for( auto& I : is2 )
        if( !hasIndex(is1, I) )
            inds.push_back(std::move(I));
    return IndexSet(inds);
    }

IndexSet
unionInds(std::vector<IndexSet> const& iss)
    {
    auto inds = std::vector<Index>();
    for( auto& is : iss )
        for(auto& I : is)
            inds.push_back(std::move(I));
    return IndexSet(inds);
    }

IndexSet
unionInds(Index const& i,
          IndexSet const& is)
    {
    auto inds = std::vector<Index>();
    inds.push_back(std::move(i));
    for( auto& I : is )
        if( I != i )
            inds.push_back(std::move(I));
    return IndexSet(inds);
    }

IndexSet
unionInds(IndexSet const& is,
          Index const& i)
    {
    auto inds = std::vector<Index>();
    for( auto& I : is )
        if( I != i )
            inds.push_back(std::move(I));
    inds.push_back(std::move(i));
    return IndexSet(inds);
    }

IndexSet
uniqueInds(IndexSet const& is1,
           IndexSet const& is2)
    {
    auto inds = std::vector<Index>();
    for( auto& I : is1 )
        if( !hasIndex(is2,I) )
            inds.push_back(std::move(I));
    return IndexSet(inds);
    }

IndexSet
uniqueInds(IndexSet const& is1,
           std::vector<IndexSet> const& is2)
    {
    return uniqueInds(is1,unionInds(is2));
    }

IndexSet
noncommonInds(IndexSet const& is1,
              IndexSet const& is2)
    {
    return unionInds(uniqueInds(is1,is2),uniqueInds(is2,is1));
    }

detail::IndexValIter
iterInds(IndexSet const& is)
    {
    return detail::IndexValIter(is);
    }

QN
flux(std::vector<IndexVal> const& ivs)
    {
    QN elt_flux;
    for(auto const& iv : ivs)
        {
        elt_flux += dir(iv)*qn(iv);
        }
    return elt_flux;
    }

#ifdef ITENSOR_USE_HDF5

void
h5_write(h5::group parent, std::string const& name, IndexSet const& is)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type","IndexSet",true);
    h5_write_attribute(g,"version",long(1));
    auto N = is.length();
    h5_write(g,"length",N);
    for(auto n : range1(N))
        {
        auto iname = format("index_%d",n);
        h5_write(g,iname,is(n));
        }
    }

void
h5_read(h5::group parent, std::string const& name, IndexSet & is)
    {
    auto g = parent.open_group(name);
    auto type = h5_read_attribute<std::string>(g,"type");
    if(type != "IndexSet") Error("Group does not contain IndexSet data in HDF5 file");
    auto N = h5_read<long>(g,"length");
    auto iv = std::vector<Index>(N);
    auto inds = IndexSetBuilder(N);
    for(auto n : range1(N))
        {
        auto iname = format("index_%d",n);
        auto i = h5_read<Index>(g,iname);
        inds.nextIndex(i);
        }
    is = inds.build();
    }

#endif

} //namespace itensor
