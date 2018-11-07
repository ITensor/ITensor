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

//
// IndexSetT Primelevel Methods
//

namespace detail {
    struct MatchInc
        {
        bool match = false;
        int inc = 0;
        MatchInc(bool m, int i) : match(m), inc(i) { }
        explicit operator bool() const { return match; };
        };
    int inline
    computeIncLast(int inc = 1) { return inc; }
    template<typename T>
    int 
    computeIncLast(const T& t, int inc) { return inc; }
    template<typename T, typename... Rest>
    int 
    computeIncLast(const T& t, Rest&&... rest) { return computeIncLast(std::forward<Rest>(rest)...); }

    template<typename Comp, typename IndexT, typename OtherT>
    MatchInc 
    computeMatchInc(const Comp& cmp,
                    const IndexT& I,
                    const OtherT& J,
                    int inc = 1)
        {
        return MatchInc(cmp(I,J),inc);
        }
    template<typename Comp, typename IndexT, typename OtherT, typename... Rest>
    MatchInc 
    computeMatchInc(const Comp& cmp,
                    const IndexT& I,
                    const OtherT& J,
                    Rest&&... rest)
        {
        if(cmp(I,J))
            {
            auto inc = computeIncLast(std::forward<Rest>(rest)...);
            return MatchInc(true,inc);
            }
        return computeMatchInc(cmp,I,std::forward<Rest>(rest)...);
        }

    template<typename Comp, typename IndexT, typename OtherT>
    bool 
    findMatch(Comp const& cmp,
              IndexT const& I,
              OtherT const& J)
        {
        return cmp(I,J);
        }
    template<typename Comp, typename IndexT, typename OtherT, typename... Rest>
    bool 
    findMatch(Comp const& cmp,
              IndexT const& I,
              OtherT const& J,
              Rest&&... rest)
        {
        if(cmp(I,J)) return true;
        return findMatch(cmp,I,std::forward<Rest>(rest)...);
        }

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

    template<typename T>
    auto
    doCheck(stdx::choice<1>,
            IndexSet const& is,
            T const& I)
        -> stdx::if_compiles_return<void,decltype(is.front()==I)>
        {
        for(auto& J : is) if(I == J) return;
        throw ITError(format("Missing index in index set\nindex = \n%s\nindex set = \n%s",I,is));
        }
    template<typename T>
    void
    doCheck(stdx::choice<2>,
            IndexSet const& is,
            T const& I)
        { }


    template<typename T>
    void
    checkHasInds(IndexSet const& is,
                 T const& I,
                 int inc = 1)
        {
        doCheck(stdx::select_overload{},is,I);
        }
    template<typename IndexT, typename T1, typename T2, typename... Rest>
    void
    checkHasInds(IndexSet const& is,
                 T1 const& I1,
                 T2 const& I2,
                 Rest&&... rest)
        {
        doCheck(stdx::select_overload{},is,I1);
        checkHasInds(is,I2,std::forward<Rest>(rest)...);
        }
} //namespace detail

//void inline
//prime(IndexSet& is, 
//      IndexType type,
//      int inc)
//    {
//    for(auto& J : is) J.prime(type,inc);
//    }

namespace detail {

struct IndexCmp
    {
    bool
    operator()(Index const& i, IndexType t) const
        {
        return i.type() == t;
        }
    bool
    operator()(IndexType t, Index const& i) const
        {
        return i.type() == t;
        }
    template<typename T1, typename T2>
    auto
    operator()(T1 const& i1, T2 const& i2) const
        -> stdx::if_compiles_return<bool,decltype(i1==i2)>
        {
        return i1 == i2;
        }
    };

} //namespace detail

template<typename... VArgs>
void 
prime(IndexSet& is, 
      VArgs&&... vargs)
    {
#ifdef DEBUG
    detail::checkHasInds(is,std::forward<VArgs>(vargs)...);
#endif
    for(auto& J : is)
        {
        auto match = 
            detail::computeMatchInc(detail::IndexCmp(),J,
                                    std::forward<VArgs>(vargs)...);
        if(match) J.prime(match.inc);
        }
#ifdef DEBUG
    detail::check(is);
#endif
    }

void inline
prime(IndexSet& is, 
      int inc) 
    { 
    prime(is,All,inc); 
    }

template<typename... VArgs>
void
primeLevel(IndexSet& is,
           VArgs&&... vargs)
    {
    constexpr size_t size = sizeof...(vargs);
    auto ints = std::array<int,size>{{static_cast<int>(vargs)...}};

    if(size != size_t(is.r()))
        {
        println("---------------------------------------------");
        println("IndexSet indices = \n",is,"\n");
        println("---------------------------------------------");
        println("Prime levels provided = ");
        for(auto& i : ints) println(i);
        println("---------------------------------------------");
        Error(format("Wrong number of prime levels passed to primeLevel (expected %d, got %d)",is.r(),size));
        }

    int i = 0;
    for(auto& J : is)
        {
        J.primeLevel(ints[i]);
        i++;
        }
    }

//template<typename IndexT, typename... Inds>
//void 
//prime(IndexSet& is, 
//      IndexT const& I1, 
//      Inds&&... rest)
//    {
//#ifdef DEBUG
//    detail::checkHasInds(is,I1,std::forward<Inds>(rest)...);
//#endif
//    auto cmp = [](const IndexT& I1, const IndexT& I2) { return I1==I2; };
//    for(auto& J : is)
//        {
//        auto match = detail::computeMatchInc(cmp,J,I1,std::forward<Inds>(rest)...);
//        if(match) J.prime(match.inc);
//        }
//#ifdef DEBUG
//    detail::check(is);
//#endif
//    }

//// Experimental prime function that
//// takes IndexVals, where the value
//// is interpreted as an increment
//template<typename IndexT, typename... IVals>
//void 
//prime(IndexSet& is,
//      const typename IndexT::indexval_type& iv1,
//      IVals&&... rest)
//    {
//#ifdef DEBUG
//    detail::checkHasInds(is,iv1,std::forward<IVals>(rest)...);
//#endif
//    int inc = 0;
//    auto cmp = [&inc](const IndexT& J, 
//                      const typename IndexT::indexval_type& iv) 
//               { 
//               inc = iv.val; 
//               return J==iv.index; 
//               };
//    for(auto& J : is)
//        {
//        auto found = detail::findMatch(cmp,J,iv1,std::forward<IVals>(rest)...);
//        if(found) J.prime(inc);
//        }
//#ifdef DEBUG
//    detail::check(is);
//#endif
//    }

template<typename... Inds>
void 
primeExcept(IndexSet& is, 
            Index const& I1, 
            Inds&&... rest)
    {
#ifdef DEBUG
    detail::checkHasInds(is,I1,std::forward<Inds>(rest)...);
#endif
    auto cmp = [](const Index& J, const Index& I) { return J==I; };
    for(auto& J : is)
        {
        auto match = detail::computeMatchInc(cmp,J,I1,std::forward<Inds>(rest)...);
        if(!match) J.prime(match.inc);
        }
#ifdef DEBUG
    detail::check(is);
#endif
    }

template<typename... ITs>
void 
primeExcept(IndexSet& is, 
            IndexType it1,
            ITs&&... rest)
    {
    auto cmp = [](const Index& J, IndexType t) { return J.type()==t; };
    for(auto& J : is)
        {
        auto match = detail::computeMatchInc(cmp,J,it1,std::forward<ITs>(rest)...);
        if(!match) J.prime(match.inc);
        }
#ifdef DEBUG
    detail::check(is);
#endif
    }

void inline
noprime(IndexSet& is, 
        IndexType type)
    {
    for(auto& J : is) J.noprime(type);
#ifdef DEBUG
    detail::check(is);
#endif
    }

template<typename... ITs>
void 
noprime(IndexSet& is,
        IndexType it1,
        IndexType it2,
        ITs&&... rest)
    {
    auto cmp = [](const Index& J, IndexType t) { return J.type()==t; };
    for(auto& J : is)
        {
        auto found = detail::findMatch(cmp,J,it1,it2,std::forward<ITs>(rest)...);
        if(found) J.noprime();
        }
#ifdef DEBUG
    detail::check(is);
#endif
	}

template<typename... Inds>
void 
noprime(IndexSet& is, 
        Index const& I1, 
        Inds&&... inds)
    {
#ifdef DEBUG
    detail::checkHasInds(is,I1,std::forward<Inds>(inds)...);
#endif
    auto cmp = [](const Index& I1, const Index& I2) { return I1==I2; };
    for(auto& J : is)
        {
        auto found = detail::findMatch(cmp,J,I1,std::forward<Inds>(inds)...);
        if(found) J.noprime();
        }
#ifdef DEBUG
    detail::check(is);
#endif
    }

void inline
mapprime(IndexSet& is, 
         int plevold, 
         int plevnew, 
         IndexType type)
	{
    for(auto& J : is) 
        {
        J.mapprime(plevold,plevnew,type);
        }
#ifdef DEBUG
    detail::check(is);
#endif
	}

namespace detail {

    bool inline
    mp_matches(Index const& I1,
               Index const& I2)
        {
        return I1.noprimeEquals(I2);
        }

    bool inline
    mp_matches(Index const& I1,
               IndexType const& T)
        {
        return I1.type() == T;
        }

    template<typename IndexOrIndexType,
             typename... VArgs>
    void
    mapprime_impl(Index & I,
                  IndexOrIndexType const& J,
                  int plev1,
                  int plev2)
        {
        if(mp_matches(I,J) && I.primeLevel()==plev1)
            {
            I.primeLevel(plev2);
            return;
            }
        }

    template<typename IndexOrIndexType,
             typename... VArgs>
    void
    mapprime_impl(Index & I,
                  IndexOrIndexType const& J,
                  int plev1,
                  int plev2,
                  VArgs&&... vargs)
        {
        if(mp_matches(I,J) && I.primeLevel()==plev1)
            {
            I.primeLevel(plev2);
            return;
            }
        mapprime_impl(I,vargs...);
        }

} //namespace detail

template<typename... VArgs>
void 
mapprime(IndexSet& is, 
         VArgs&&... vargs)
	{
    static_assert(sizeof...(vargs)%3==0,
                  "Wrong number of arguments to mapprime");
    for(auto& I : is)
        {
        detail::mapprime_impl(I,vargs...);
        }
#ifdef DEBUG
    detail::check(is);
#endif
	}

void inline
sim(IndexSet & is, 
    IndexType t)
    {
    for(auto n : range(is))
        {
        if(is[n].type() == t)
            {
            is[n] = sim(is[n]);
            }
        }
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

//
//
// IndexSet helper methods
//
//


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
findindex(const IndexSet& iset, 
          const Index& I)
    {
    for(long j = 0; j < iset.r(); ++j)
        {
        if(iset[j] == I) return j;
        }
    return -1;
    }

Index inline
findtype(IndexSet const& iset, IndexType t)
	{
    for(auto& J : iset)
        {
        if(J.type() == t) return J;
        }
    Error("findtype failed."); 
    return Index();
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
hasindex(const IndexSet& iset, 
         const Index& I)
	{
    for(long j = 0; j < iset.r(); ++j)
        {
        if(iset[j] == I) return true;
        }
    return false;
	}

bool inline
hastype(const IndexSet& iset, 
        IndexType t)
	{
    for(const auto& J : iset)
        {
        if(J.type() == t) return true;
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
            newind.nextIndStr(Lis[j],Lind[j]);
        }
    for(decltype(Ris.r()) j = 0; j < Ris.r(); ++j)
        {
        if(Rind[j] > 0) //uncontracted
            newind.nextIndStr(Ris[j],Rind[j]);
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

//inline std::ostream&
//operator<<(std::ostream& s, IndexSet const& is)
//    {
//    auto size = is.size();
//    if(size > 0) s << is[0];
//    for(decltype(size) j = 1; j < size; ++j)
//        {
//        s << "\n" << is[j];
//        }
//    return s;
//    }

inline std::ostream&
operator<<(std::ostream& s, IndexSet const& is)
    {
    for(auto i : range1(is.r()))
        { 
        s << is.index(i);
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

} //namespace itensor

#endif
