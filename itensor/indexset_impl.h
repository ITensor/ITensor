//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_IMPL_H
#define __ITENSOR_INDEXSET_IMPL_H

namespace itensor {

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


void inline IndexSet::
setPrime(int plnew, TagSet const& tsmatch) 
    { 
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.setPrime(plnew); 
    }

void inline IndexSet::
setPrime(int plnew, Index const& imatch) 
    { 
    for(auto& J : *this) if(J==imatch) J.setPrime(plnew); 
    }

void inline IndexSet::
noPrime(TagSet const& tsmatch) 
    { 
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.setPrime(0); 
    }

void inline IndexSet::
noPrime(Index const& imatch) 
    { 
    for(auto& J : *this) if(J==imatch) J.setPrime(0); 
    }

// TODO: should this just be replaced by setPrime()?
void inline IndexSet::
mapPrime(int plold, 
         int plnew, 
         TagSet const& tsmatch) 
    { 
    for(auto& J : *this) if(matchTagsPrime(J,tsmatch,plold)) J.setPrime(plnew); 
    }

void inline IndexSet::
prime(int plinc, 
      TagSet const& tsmatch) 
    { 
    for(auto& J : *this) if(hasTags(J,tsmatch)) J.prime(plinc);
    }

void inline IndexSet::
prime(TagSet const& tsmatch) 
    { 
    this->prime(1,tsmatch);
    }

void inline IndexSet::
prime(int plinc, 
      Index const& imatch) 
    { 
    for(auto& J : *this) if(J==imatch) J.prime(plinc);
    }

void inline IndexSet::
prime(Index const& imatch) 
    { 
    this->prime(1,imatch);
    }

void inline IndexSet::
swapPrime(int pl1, 
          int pl2, 
          TagSet const& tsmatch)
    {
    int tempLevel = 99999;
#ifdef DEBUG
    for(auto& I : *this)
        {
        if(I.primeLevel() == tempLevel)
            {
            println("tempLevel = ",tempLevel);
            Error("swapPrime fails if an index has primeLevel==tempLevel");
            }
        }
#endif
    this->mapPrime(pl1,tempLevel,tsmatch);
    this->mapPrime(pl2,pl1,tsmatch);
    this->mapPrime(tempLevel,pl2,tsmatch);
    }

void inline IndexSet::
replaceTags(TagSet const& tsold, 
            TagSet const& tsnew, 
            TagSet const& tsmatch, 
            int plmatch)
    {
    for(auto& J : *this)
        {
        if(matchTagsPrime(J,tsmatch,plmatch) && hasTags(J,tsold))
            {
            J.removeTags(tsold);
            J.addTags(tsnew);
            }
        }
    }

void inline IndexSet::
swapTags(TagSet const& ts1, 
         TagSet const& ts2, 
         TagSet const& tsmatch, 
         int plmatch)
    {
    auto tmp = TagSet("df4sd321");  // Some random temporary tag. TODO: Debug level check that ts1!=tmp && ts2!=tmp
    this->replaceTags(ts1,tmp,tsmatch,plmatch);
    this->replaceTags(ts2,ts1,tsmatch,plmatch);
    this->replaceTags(tmp,ts2,tsmatch,plmatch);
    }

void inline IndexSet::
setTags(TagSet const& tsnew, 
        TagSet const& tsmatch, 
        int plmatch) 
    { 
    for(auto& J : *this) if(matchTagsPrime(J,tsmatch,plmatch)) J.setTags(tsnew); 
    }

void inline IndexSet::
addTags(TagSet const& tsadd, 
        TagSet const& tsmatch, int plmatch) 
    { 
    for(auto& J : *this) if(matchTagsPrime(J,tsmatch,plmatch)) J.addTags(tsadd); 
    }

void inline IndexSet::
addTags(TagSet const& tsadd, 
        int plmatch) 
    { 
    for(auto& J : *this) if(J.primeLevel()==plmatch) J.addTags(tsadd); 
    }


void inline IndexSet::
addTags(TagSet const& tsadd, 
        Index const& imatch) 
    { 
    for(auto& J : *this) if(J==imatch) J.addTags(tsadd);
    }

void inline IndexSet::
removeTags(TagSet const& tsremove, 
           TagSet const& tsmatch,
           int plmatch) 
    { 
    for(auto& J : *this) if(matchTagsPrime(J,tsmatch,plmatch)) J.removeTags(tsremove); 
    }

void inline IndexSet::
removeTags(TagSet const& tsremove, 
           int plmatch) 
    { 
    for(auto& J : *this) if(J.primeLevel()==plmatch) J.removeTags(tsremove);
    }

void inline IndexSet::
removeTags(TagSet const& tsremove, 
           Index const& imatch) 
    { 
    for(auto& J : *this) if(J==imatch) J.removeTags(tsremove); 
    }

void inline IndexSet::
removeQNs()
    {
    for(auto& J : *this) J.removeQNs();
    }

//
// Methods for Manipulating IndexSetT
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

namespace detail {
    void inline
    check(IndexSet const& is)
        {
        //Check if any duplicate indices
        for(size_t j = 0; j < is.size(); ++j) 
        for(size_t k = 0; k < is.size(); ++k)
            if(k != j && is[j] == is[k])
                {
                println("index set = \n",is);
                throw ITError("Duplicate indices in index set");
                }
        }

} //namespace detail

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
    Error("dir: Index not found");
    return In;
    }


Index inline
finddir(IndexSet const& iset, Arrow dir)
    {
    for(const auto& J : iset)
        {
        if(J.dir() == dir) return J;
        }
    Error("Couldn't find index with specified dir");
    return Index();
    }

//
// Given IndexSet iset and Index I,
// return int j such that iset[j] == I.
// If not found, returns -1
//
long inline
findIndex(const IndexSet& iset, 
          const Index& I)
    {
    for(long j = 0; j < iset.r(); ++j)
        {
        if(iset[j] == I) return j;
        }
    return -1;
    }


////
//// Compute the permutation P taking an IndexSetT iset
//// to oset (of type IndexSetT or array<Index,NMAX>)
////
//template <class IndexT>
//void
//getperm(const IndexSet& iset, 
//        const typename IndexSet::storage& oset, 
//        Permutation& P)
//	{
//	for(int j = 0; j < iset.r(); ++j)
//	    {
//	    bool got_one = false;
//	    for(int k = 0; k < iset.r(); ++k)
//            {
//            if(oset[j] == iset[k])
//                { 
//                P.setFromTo(j+1,k+1); 
//                got_one = true; 
//                break;
//                }
//            }
//	    if(!got_one)
//            {
//            println("j = ",j);
//            println("iset =");
//            for(int j = 0; j < iset.r(); ++j)
//                printfln("%d %s",j,iset[j]);
//            println("\noset = ");
//            for(int j = 0; j < iset.r(); ++j)
//                printfln("%d %s",j,oset[j]);
//            println();
//            //printfln("iset uniqueReal = %.15E",iset.uniqueReal());
//            //Real our = 0;
//            //for(int i = 0; i < iset.r(); ++i)
//            //    {
//            //    our += oset[i].uniqueReal();
//            //    }
//            //printfln("oset uniqueReal = %.15E",our);
//            //printfln("uniqueReal diff = %.15E",fabs(our-iset.uniqueReal()));
//            throw ITError("IndexSetT::getperm: no matching index");
//            }
//	    }
//	}

bool inline
hasIndex(const IndexSet& iset, 
         const Index& I)
	{
    for(long j = 0; j < iset.r(); ++j)
        {
        if(iset[j] == I) return true;
        }
    return false;
	}

long inline
minM(const IndexSet& iset)
    {
    if(iset.empty()) return 1l;
    auto mm = iset[0].m();
    for(long j = 1; j < iset.r(); ++j)
        mm = std::min(mm,iset[j].m());

    return mm;
    }

long inline
maxM(const IndexSet& iset)
    {
    if(iset.empty()) return 1l;

    auto mm = iset[0].m();
    for(long j = 1; j < iset.r(); ++j)
        mm = std::max(mm,iset[j].m());

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
    auto nuniq = Lis.r()+Ris.r()-2*ncont;
    auto newind = RangeBuilderT<IndexSet>(nuniq);
    //Below we are "cheating" and using the .str
    //field of each member of newind to hold the 
    //labels which will go into Nind so they will
    //be sorted along with the .ext members (the indices of Nis)
    //Later we will set them back to zero
    //IndexSetT constructor anyway
    for(decltype(Lis.r()) j = 0; j < Lis.r(); ++j)
        {
        if(Lind[j] > 0) //uncontracted
            {
            if(LremoveQNs) newind.nextIndStr(removeQNs(Lis[j]),Lind[j]);
            else newind.nextIndStr(Lis[j],Lind[j]);
            }
        }
    for(decltype(Ris.r()) j = 0; j < Ris.r(); ++j)
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
    computeLabels(Lis,Lis.r(),Ris,Ris.r(),Lind,Rind);
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
    auto nuniq = Lis.r()+Ris.r()-2*nmerge;
    auto NisBuild = RangeBuilderT<IndexSet>(nmerge+nuniq);
    Nind.resize(nmerge+nuniq);
    long n = 0;

    for(decltype(Lis.r()) j = 0; j < Lis.r(); ++j)
        {
        if(Lind[j] < 0) //merged
            {
            NisBuild.nextIndex(Lis[j]);
            Nind[n++] = Lind[j];
            }
        }
    for(decltype(Lis.r()) j = 0; j < Lis.r(); ++j)
        {
        if(Lind[j] > 0) //unmerged/unique
            {
            NisBuild.nextIndex(Lis[j]);
            Nind[n++] = Lind[j];
            }
        }
    for(decltype(Ris.r()) j = 0; j < Ris.r(); ++j)
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
    for(auto i : range1(is.r()))
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
    if(is.r() > 0 && hasQNs(is.front()))
        {
        for(long n = 1; n < is.r(); n += 1)
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
          TagSet const& tsmatch, 
          int plmatch)
    {
    for(auto& J : is) if(matchTagsPrime(J,tsmatch,plmatch)) return J;
    //TODO: make this a debug error
    Error("No index with those tags and prime level found");
    return Index();
    }

Index inline
findIndexExact(IndexSet const& is,
               TagSet const& tsmatch, 
               int plmatch)
    {
    auto j = Index();
    for(auto& J : is)
        {
        if(matchTagsPrimeExact(J,tsmatch,plmatch))
            {
            if(j) Error("Multiple indices with those tags and prime level found");
            else j = J;
            }
        }
    if(!j) Error("No index with those tags and prime level found");
    return j;
    }


} //namespace itensor

#endif
