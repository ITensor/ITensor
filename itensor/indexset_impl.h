//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_IMPL_H
#define __ITENSOR_INDEXSET_IMPL_H

namespace itensor {

namespace detail {

void inline
check(IndexSet const& is)
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

void inline
checkIndexPositions(std::vector<int> const& is)
    {
    //Check if there any duplicate index positions
    //(it is ok if Indices are not found, they are just ignored)
    for(size_t j = 0; j < is.size(); ++j) 
    for(size_t k = j+1; k < is.size(); ++k)
        if( (is[j] != -1) && (is[j] == is[k]) )
            throw ITError("An Index was found more than once in the IndexSet");
    }

} //namespace detail

void inline IndexSet::
dag() { for(auto& J : *this) J.dag(); }

IndexSet::iterator inline IndexSet::
begin() { return iterator{*this}; }

IndexSet::iterator inline IndexSet::
end() { return iterator::makeEnd(*this); }

IndexSet::const_iterator inline IndexSet::
begin() const { return const_iterator{*this}; }

IndexSet::const_iterator inline IndexSet::
end() const { return const_iterator::makeEnd(*this); }

IndexSet::const_iterator inline IndexSet::
cbegin() const { return begin(); }

IndexSet::const_iterator inline IndexSet::
cend() const { return end(); }

//
// Tag functions
//

void inline IndexSet::
setTags(TagSet const& tsnew)
    {
    for(auto& J : *this) J.setTags(tsnew);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
setTags(TagSet const& tsnew, 
        IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).setTags(tsnew);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
setTags(TagSet const& tsnew, 
        TagSet const& tsmatch)
    { 
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.setTags(tsnew); 
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
addTags(TagSet const& tsadd)
    {
    for(auto& J : *this) J.addTags(tsadd);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
addTags(TagSet const& tsadd, 
        IndexSet const& ismatch) 
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).addTags(tsadd);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
addTags(TagSet const& tsadd, 
        TagSet const& tsmatch)
    { 
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.addTags(tsadd); 
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
removeTags(TagSet const& tsremove)
    {
    for(auto& J : *this) J.removeTags(tsremove);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
removeTags(TagSet const& tsremove, 
           IndexSet const& ismatch) 
    { 
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).removeTags(tsremove);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
removeTags(TagSet const& tsremove, 
           TagSet const& tsmatch)
    { 
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.removeTags(tsremove); 
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew)
    {
    for(auto& J : *this) J.replaceTags(tsold,tsnew);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
replaceTags(TagSet const& tsold, 
            TagSet const& tsnew, 
            IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).replaceTags(tsold,tsnew);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
replaceTags(TagSet const& tsold,
            TagSet const& tsnew,
            TagSet const& tsmatch)
    {
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.replaceTags(tsold,tsnew);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

template<typename... VarArgs>
void IndexSet::
swapTags(TagSet const& ts1, 
         TagSet const& ts2, 
         VarArgs&&... vargs)
    {
    TagSet tempTags;
    if(primeLevel(ts1) < 0)
      {
      if(primeLevel(ts2) >= 0) Error("Error in swapTags: cannot swap integer tag with non-integer tag");
      tempTags = TagSet("df4sd32");
      }
    else
      {
      if(primeLevel(ts2) < 0) Error("Error in swapTags: cannot swap integer tag with non-integer tag");
      tempTags = TagSet("df4sd32,431543");
      }
#ifdef DEBUG
    for(auto& I : *this)
        {
        if(hasTags(I,tempTags))
            {
            println("tempTags = ",tempTags);
            throw ITError("swapTags fails if an index has tags tempTags");
            }
        }
#endif
    this->replaceTags(ts1,tempTags,std::forward<VarArgs>(vargs)...);
    this->replaceTags(ts2,ts1,std::forward<VarArgs>(vargs)...);
    this->replaceTags(tempTags,ts2);
    }

//
// Integer tag convenience functions
//

void inline IndexSet::
setPrime(int plnew)
    { 
    for(auto& J : *this) J.setPrime(plnew); 
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
setPrime(int plnew, IndexSet const& ismatch)
    {
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).setPrime(plnew);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
setPrime(int plnew, TagSet const& tsmatch)
    { 
    for(auto& J : *this) if( hasTags(J,tsmatch) ) J.setPrime(plnew); 
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
prime(int plinc)
    { 
    for(auto& J : *this) J.prime(plinc);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

template<typename... VarArgs>
void IndexSet::
prime(int plinc, IndexSet const& ismatch)
    { 
    auto ilocs = indexPositions(*this,ismatch);
    for(auto i : ilocs) parent::index(i).prime(plinc);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
prime(int plinc, TagSet const& tsmatch)
    { 
    for(auto& J : *this) if( hasTags(J,tsmatch) ) J.prime(plinc);
#ifdef DEBUG
    detail::check(*this);
#endif
    }

void inline IndexSet::
removeQNs()
    {
    for(auto& J : *this) J.removeQNs();
    }

//
// Methods for Manipulating IndexSet
//

void inline
write(std::ostream& s, IndexSet const& is)
    {
    using parent = typename IndexSet::parent;
    parent const& pr = is;
    itensor::write(s,pr);
    }

void inline
read(std::istream& s, IndexSet & is)
    {
    using parent = typename IndexSet::parent;
    parent & pr = is;
    itensor::read(s,pr);
    }

void inline
sim(IndexSet & is, 
    Index const& I)
    {
    for(auto n : range(is))
        {
        if(is[n] == I)
            {
            is[n] = sim(I);
            }
        }
    }

Arrow inline
dir(const IndexSet& is, const Index& I)
    {
    for(const auto& J : is)
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
int inline
indexPosition(IndexSet const& is, 
              Index const& imatch)
    {
    for(auto j : range(order(is))) 
        if( is[j] == imatch ) return j;
    return -1;
    }


std::vector<int> inline
indexPositions(IndexSet const& is,
               IndexSet const& ismatch)
    {
    auto ilocs = std::vector<int>();
    for(const auto& J : ismatch)
        {
        auto loc = indexPosition(is,J);
        if( loc != -1 ) ilocs.push_back(loc);
        }
#ifdef DEBUG
    detail::checkIndexPositions(ilocs);
#endif
    return ilocs;
    }

bool inline
hasIndex(const IndexSet& iset, 
         const Index& I)
	{
  for(long j = 0; j < iset.order(); ++j)
      if(iset[j] == I) return true;
  return false;
	}

long inline
minDim(const IndexSet& iset)
    {
    if(iset.empty()) return 1l;
    auto mm = dim(iset[0]);
    for(long j = 1; j < iset.order(); ++j)
        mm = std::min(mm,dim(iset[j]));

    return mm;
    }

long inline
maxDim(const IndexSet& iset)
    {
    if(iset.empty()) return 1l;

    auto mm = dim(iset[0]);
    for(long j = 1; j < iset.order(); ++j)
        mm = std::max(mm,dim(iset[j]));

    return mm;
    }

template<class LabelT>
void
contractIS(IndexSet const& Lis,
           LabelT const& Lind,
           IndexSet const& Ris,
           LabelT const& Rind,
           IndexSet & Nis,
           LabelT & Nind,
           bool sortResult)
    {
    //Don't allow mixing of Indices with and without QNs
    auto LhasQNs = hasQNs(Lis);
    auto RhasQNs = hasQNs(Ris);
    auto LremoveQNs = false;
    auto RremoveQNs = false;
    if(LhasQNs && !RhasQNs) LremoveQNs = true;
    else if(!LhasQNs && RhasQNs) RremoveQNs = true;

    long ncont = 0;
    for(auto& i : Lind) if(i < 0) ++ncont;
    auto nuniq = Lis.order()+Ris.order()-2*ncont;
    auto newind = RangeBuilderT<IndexSet>(nuniq);
    //Below we are "cheating" and using the .str
    //field of each member of newind to hold the 
    //labels which will go into Nind so they will
    //be sorted along with the .ext members (the indices of Nis)
    //Later we will set them back to zero
    //IndexSetT constructor anyway
    for(decltype(Lis.order()) j = 0; j < Lis.order(); ++j)
        {
        if(Lind[j] > 0) //uncontracted
            {
            if(LremoveQNs) newind.nextIndStr(removeQNs(Lis[j]),Lind[j]);
            else newind.nextIndStr(Lis[j],Lind[j]);
            }
        }
    for(decltype(Ris.order()) j = 0; j < Ris.order(); ++j)
        {
        if(Rind[j] > 0) //uncontracted
            {
            if(RremoveQNs) newind.nextIndStr(removeQNs(Ris[j]),Rind[j]);
            else newind.nextIndStr(Ris[j],Rind[j]);
            }
        }
    if(sortResult) newind.sortByIndex();
    Nind.resize(newind.size());
    for(decltype(newind.size()) j = 0; j < newind.size(); ++j)
        {
        Nind[j] = newind.stride(j);
        }
    Nis = newind.build();
    Nis.computeStrides();
    }

void inline
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

template<class LabelT>
void
ncprod(IndexSet const& Lis,
       LabelT const& Lind,
       IndexSet const& Ris,
       LabelT const& Rind,
       IndexSet & Nis,
       LabelT & Nind)
    {
    long nmerge = 0;
    for(auto& i : Lind) if(i < 0) ++nmerge;
    auto nuniq = Lis.order()+Ris.order()-2*nmerge;
    auto NisBuild = RangeBuilderT<IndexSet>(nmerge+nuniq);
    Nind.resize(nmerge+nuniq);
    long n = 0;

    for(decltype(Lis.order()) j = 0; j < Lis.order(); ++j)
        {
        if(Lind[j] < 0) //merged
            {
            NisBuild.nextIndex(Lis[j]);
            Nind[n++] = Lind[j];
            }
        }
    for(decltype(Lis.order()) j = 0; j < Lis.order(); ++j)
        {
        if(Lind[j] > 0) //unmerged/unique
            {
            NisBuild.nextIndex(Lis[j]);
            Nind[n++] = Lind[j];
            }
        }
    for(decltype(Ris.order()) j = 0; j < Ris.order(); ++j)
        {
        if(Rind[j] > 0) //unmerged/unique
            {
            NisBuild.nextIndex(Ris[j]);
            Nind[n++] = Rind[j];
            }
        }
    Nis = NisBuild.build();
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

bool inline
hasQNs(IndexSet const& is)
    {
    for(auto& I : is) if(hasQNs(I)) return true;
    return false;
    }

bool inline
hasQNs(std::vector<Index> const& inds)
    {
    for(auto& I : inds) if(hasQNs(I)) return true;
    return false;
    }

void inline
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

Index inline
findIndex(IndexSet const& is,
          TagSet const& tsmatch)
    {
    for(auto& J : is) if(hasTags(J,tsmatch))
        return J;
    return Index();
    }

} //namespace itensor

#endif
