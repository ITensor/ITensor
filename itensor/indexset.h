//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include <algorithm>
#include "itensor/util/readwrite.h"
#include "itensor/util/safe_ptr.h"
#include "itensor/tensor/range.h"
#include "itensor/tensor/permutation.h"
#include "itensor/index.h"


namespace itensor {

template <class IndexT>
class IndexSetT;

class IQIndex;

//
// IndexSetT
//
// Aliases:
using IndexSet = IndexSetT<Index>;
using IQIndexSet = IndexSetT<IQIndex>;

//
// When constructed from a collection of indices,
// (as an explicit set of arguments or via
// a container) puts the indices with m>1 to the
// front and those with m==1 at the back but otherwise
// keeps the indices in the order given.
//

template <class IndexT>
class IndexSetT
    {
    public:

    using storage = std::vector<IndexT>;
    using value_type = IndexT;
    using iterator = typename storage::iterator;
    using const_iterator = typename storage::const_iterator;
    using index_type = IndexT;
    using indexval_type = typename IndexT::indexval_type;

    IndexSetT();

    explicit
    IndexSetT(const IndexT& i1);

    // construct from 2 or more indices
    template <typename... Inds>
    IndexSetT(const IndexT& i1, 
              const IndexT& i2,
              const Inds&... inds)
        {
        init(std::array<IndexT,2+sizeof...(inds)>{{i1,i2,inds...}});
        }

    explicit
    IndexSetT(const std::vector<IndexT>& ii) { init(ii); }

    template<size_t N>
    explicit
    IndexSetT(const std::array<IndexT,N>& ii) { init(ii); }

    IndexSetT(std::initializer_list<IndexT> ii) { init(ii); }
    
    explicit operator bool() const { return !index_.empty(); }

    //
    // Accessor Methods
    //

    long
    dim(long i) const;

    long
    stride(long i) const;

    int
    r() const { return index_.size(); }
    
    size_t
    size() const { return index_.size(); }

    bool
    empty() const { return index_.empty(); }

    int
    rn() const { return rn_; }

    // 0-indexed access
    const IndexT&
    operator[](int j) const;

    IndexT&
    operator[](int j);

    // 1-indexed access
    const IndexT&
    index(int j) const;

    const IndexT&
    front() const { return index_.front(); }

    const IndexT&
    back() const { return index_.back(); }

    const_iterator
    begin() const { return index_.begin(); }

    const_iterator
    end() const { return index_.end(); }

    iterator
    begin() { return index_.begin(); }

    iterator
    end() { return index_.end(); }

    const_iterator
    cbegin() const { return index_.cbegin(); }

    const_iterator
    cend() const { return index_.cend(); }


    //
    // Other Methods
    //

    void 
    addindex(const IndexT& I);

    void
    replaceIndex(const IndexT& oind, const IndexT& nind);

    void
    swap(IndexSetT& other);

    void
    clear();

    void
    dag();

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    private:

    //////////

    storage index_;
    std::vector<long> stride_;
    int rn_;

    /////////

    template<class Iterable>
    void
    init(Iterable&& inds);

    };

//
// IndexSetT Primelevel Methods
//

template<typename IndexT, typename... Types>
void 
prime(IndexSetT<IndexT>& is, 
      IndexType type,
      int inc = 1);

template<typename IndexT, typename... Types>
void 
prime(IndexSetT<IndexT>& is, 
      IndexType type1,
      Types&&... rest);

template<typename IndexT>
void 
prime(IndexSetT<IndexT>& is, int inc = 1) { prime(is,All,inc); }

template<typename IndexT, typename... IVals>
void 
prime(IndexSetT<IndexT>& is,
      const typename IndexT::indexval_type& iv1,
      IVals&&... ivs);

//
//Given a list of indices and an increment (an int)
//as the optional last argument (default is inc=1)
//increment all indices NOT listed in the arguments
//by the amout inc.
//
//For example, primeExcept(is,I1,I3,I4,I7,2);
//will increment all prime levels by 2 except for
//those of I1,I3,I4, and I7.
//
template<typename IndexT, typename... Inds>
void 
primeExcept(IndexSetT<IndexT>& is, 
            const IndexT& I1, 
            Inds&&... inds);

template<typename IndexT, typename... ITs>
void 
primeExcept(IndexSetT<IndexT>& is, 
            IndexType it,
            ITs&&... etc);

template<typename IndexT>
void 
noprime(IndexSetT<IndexT>& is, IndexType type = All);

template<typename IndexT, typename... ITs>
void 
noprime(IndexSetT<IndexT>& is,
        IndexType it1,
        IndexType it2,
        ITs&&... rest);

template<typename IndexT, typename... Inds>
void 
noprime(IndexSetT<IndexT>& is, 
        const IndexT& I1, 
        Inds&&... inds);

template<typename IndexT>
void 
mapprime(IndexSetT<IndexT>& is, int plevold, int plevnew, IndexType type = All);

//
// IndexSetT helper methods
//


template<class IndexT>
Arrow
dir(const IndexSetT<IndexT>& is, const IndexT& I);


template <class IndexT>
const IndexT&
finddir(const IndexSetT<IndexT>& iset, Arrow dir);

//
// Given IndexSetT<IndexT> iset and IndexT I,
// return int j such that iset[j] == I.
// If not found, returns -1
//
template <class IndexT>
int
findindex(const IndexSetT<IndexT>& iset, 
          const IndexT& I);

template <class IndexT>
const IndexT&
findtype(const IndexSetT<IndexT>& iset, IndexType t);

//
// Compute the permutation P taking an IndexSetT iset
// to oset (of type IndexSetT or array<IndexT,NMAX>)
//
template <class IndexT>
void
getperm(const IndexSetT<IndexT>& iset, 
        const typename IndexSetT<IndexT>::storage& oset, 
        Permutation& P);

template <class IndexT>
bool
hasindex(const IndexSetT<IndexT>& iset, 
         const IndexT& I);

template <class IndexT>
bool
hastype(const IndexSetT<IndexT>& iset, 
        IndexType t);

template <class IndexT>
long
minM(const IndexSetT<IndexT>& iset);

template <class IndexT>
long
maxM(const IndexSetT<IndexT>& iset);

template<class IndexT, class ContainerT>
void
contractIS(const IndexSetT<IndexT>& Lis,
          const ContainerT& Lind,
          const IndexSetT<IndexT>& Ris,
          const ContainerT& Rind,
          IndexSetT<IndexT>& Nis,
          bool sortResult);

template <class IndexT>
std::ostream&
operator<<(std::ostream& s, const IndexSetT<IndexT>& is);

} //namespace itensor

#include "itensor/indexset.ih"

#endif
