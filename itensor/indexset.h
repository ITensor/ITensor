//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include "index.h"
#include "permutation.h"

//
// IndexSet
//
class IndexSet
    {
    public:

    IndexSet();

    explicit
    IndexSet(const Index& i1);

    IndexSet(const Index& i1, const Index& i2);

    IndexSet(Index i1, Index i2, Index i3,
             Index i4 = Index::Null(), 
             Index i5 = Index::Null(), 
             Index i6 = Index::Null(),
             Index i7 = Index::Null(), 
             Index i8 = Index::Null());

    template <class Iterable>
    IndexSet(const Iterable& ii, int size, int& alloc_size, int offset = 0);

    IndexSet(const IndexSet& other, const Permutation& P);

    //
    // Accessor Methods
    //

    int
    r() const { return r_; }

    int
    rn() const { return rn_; }

    const Index&
    index(int j) const { return GET(index_,j); }

    int
    m(int j) const { return GET(index_,j).m(); }

    typedef boost::array<Index,NMAX+1>::const_iterator 
    index_it;

    //Can be used for iteration over Indices in a Foreach loop
    //e.g. Foreach(const Index& I, t.index() ) { ... }
    const std::pair<index_it,index_it> 
    index() const  
        { return std::make_pair(index_.begin()+1,index_.begin()+r_+1); }

    Real
    uniqueReal() const { return ur_; }

    //
    // Index Analysis
    //

    const Index&
    findtype(IndexType t) const;

    int 
    findindex(const Index& I) const;

    int 
    findindexn(const Index& I) const;

    int 
    findindex1(const Index& I) const;

    bool 
    has_common_index(const IndexSet& other) const;
    
    bool 
    hasindex(const Index& I) const;

    bool 
    hasindexn(const Index& I) const;

    bool 
    hasindex1(const Index& I) const;

    bool
    hasAllIndex(const boost::array<Index,NMAX+1>& I, int nind) const;

    void
    getperm(const boost::array<Index,NMAX+1>& I, Permutation& P) const;

    void
    getperm(const IndexSet& other, Permutation& P) const;

    int
    minM() const;

    int
    maxM() const;

    //
    // Primelevel Methods
    //

    void 
    noprime(IndexType type = All);

    void 
    doprime(IndexType type, int inc = 1);

    void 
    primeall() { doprime(All,1); }

    void 
    primesite() { doprime(Site,1); }

    void 
    primelink() { doprime(Link,1); }

    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    void 
    mapprimeind(const Index& I, int plevold, int plevnew, 
                IndexType type = All);

    void 
    primeind(const Index& I, int inc = 1)
        { mapindex(I,primed(I,inc)); }

    void 
    primeind(const Index& I, const Index& J);

    void 
    noprimeind(const Index& I) { mapindex(I,deprimed(I)); }

    //
    // Methods for Manipulating IndexSets
    //
    // Warning: these can overwrite other
    // Indices if not used properly
    //

    void 
    mapindex(const Index& i1, const Index& i2);

    void 
    addindex(const Index& I);

    void 
    addindexn(const boost::array<Index,NMAX+1>& indices, int n);

    void 
    addindexn(const Index& I);

    void 
    addindex1(const boost::array<Index,NMAX+1>& indices, int n);

    void 
    addindex1(const std::vector<Index>& indices);

    void 
    addindex1(const Index& I);

    //Removes the jth index as found by findindex
    void 
    removeindex1(int j);

    void 
    removeindex1(const Index& I) 
        { removeindex1(findindex1(I)); }

    void
    setUniqueReal();

    void
    swap(IndexSet& other);

    void
    clear();

    //
    // Other Methods
    //

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    friend std::ostream&
    operator<<(std::ostream& s, const IndexSet& is);

    //////////
    //
    // Data Members
    //

    boost::array<Index,NMAX+1> index_;

    int rn_,
        r_;

    Real ur_;

    //
    /////////

    };

template <class Iterable>
IndexSet::
IndexSet(const Iterable& ii, int size, int& alloc_size, int offset)
    :
    r_(size)
    { 
    sortIndices(ii,size,rn_,alloc_size,index_,offset);
    setUniqueReal();
    }


template<class Iterable>
void
sortIndices(const Iterable& I, int ninds, int& rn_, int& alloc_size, 
            boost::array<Index,NMAX+1>& index_, int offset = 0)
    {
    assert(ninds <= NMAX);

    rn_ = 0;
    alloc_size = 1;

    int r1_ = 0;
    boost::array<const Index*,NMAX+1> index1_;

    for(int n = offset; n < ninds+offset; ++n)
        {
        const Index& i = I[n];
        DO_IF_DEBUG(if(i == Index::Null()) Error("Null Index in sortIndices");)
        if(i.m()==1) 
            { index1_[++r1_] = &i; }
        else         
            { 
            index_[++rn_] = i; 
            alloc_size *= i.m(); 
            }
        }
    for(int l = 1; l <= r1_; ++l) 
        index_[rn_+l] = *(index1_[l]);
    }

#endif
