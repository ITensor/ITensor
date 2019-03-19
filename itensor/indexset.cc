//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
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

void IndexSet::
dag() { for(auto& J : *this) J.dag(); }

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

void IndexSet::
setTags(TagSet const& tsnew)
    {
    for(auto& J : *this) J.setTags(tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
setTags(TagSet const& tsnew,
        IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).setTags(tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
setTags(TagSet const& tsnew,
        TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.setTags(tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
addTags(TagSet const& tsadd)
    {
    for(auto& J : *this) J.addTags(tsadd);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
addTags(TagSet const& tsadd,
        IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).addTags(tsadd);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
addTags(TagSet const& tsadd,
        TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.addTags(tsadd);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
removeTags(TagSet const& tsremove)
    {
    for(auto& J : *this) J.removeTags(tsremove);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
removeTags(TagSet const& tsremove,
           IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).removeTags(tsremove);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
removeTags(TagSet const& tsremove,
           TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.removeTags(tsremove);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew)
    {
    for(auto& J : *this) J.replaceTags(tsold,tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew,
            IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).replaceTags(tsold,tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew,
            TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.replaceTags(tsold,tsnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
setPrime(int plnew)
    {
    for(auto& J : *this) J.setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
setPrime(int plnew, IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
setPrime(int plnew, TagSet const& tsmatch)
    {
    for(auto& J : *this) if( hasTags(J,tsmatch) ) J.setPrime(plnew);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
prime(int plinc)
    {
    for(auto& J : *this) J.prime(plinc);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
prime(int plinc, IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).prime(plinc);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

void IndexSet::
prime(int plinc, TagSet const& tsmatch)
    {
    for(auto& J : *this) if( hasTags(J,tsmatch) ) J.prime(plinc);
#ifdef DEBUG
    checkIndexSet(*this);
#endif
    }

IndexSet
sim(IndexSet is)
    {
    for(auto& J : is)
        J = sim(J);
    return is;
    }

void IndexSet::
removeQNs()
    {
    for(auto& J : *this) J.removeQNs();
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
hasIndex(IndexSet const& iset,
         Index const& I)
  {
  for(long j = 0; j < iset.order(); ++j)
      if(iset[j] == I) return true;
  return false;
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

inline std::ostream&
operator<<(std::ostream& s, IndexSet const& is)
    {
    for(auto i : range1(is.order()))
        {
        if(hasQNs(is)) s << is.index(i) << std::endl;
        else s << is.index(i) << " ";
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

Index
findIndex(IndexSet const& is,
          TagSet const& tsmatch)
    {
    for(auto& J : is) if(hasTags(J,tsmatch))
        return J;
    return Index();
    }

IndexSet
unionInds(IndexSet const& is1,
          IndexSet const& is2)
    {
    auto inds = std::vector<Index>();
    for( auto J : is1 )
        inds.push_back(J);
    for( auto J : is2 ) if( !hasIndex(is1, J) )
        inds.push_back(J);
    return IndexSet(inds);
    }

} //namespace itensor
