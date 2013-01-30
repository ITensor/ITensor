//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQINDEXSET_H
#define __ITENSOR_IQINDEXSET_H
#include "iqindex.h"

//
// IQIndexSet
//
class IQIndexSet
    {
    public:

    IQIndexSet();

    explicit
    IQIndexSet(const IQIndex& i1);

    IQIndexSet(const IQIndex& i1, const IQIndex& i2);

    IQIndexSet(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3);

    IQIndexSet(IQIndex i1, IQIndex i2, 
               IQIndex i3, IQIndex i4,
               IQIndex i5 = IQIndex::Null(), 
               IQIndex i6 = IQIndex::Null(),
               IQIndex i7 = IQIndex::Null(), 
               IQIndex i8 = IQIndex::Null());

    IQIndexSet(std::vector<IQIndex>& iqinds);

    IQIndexSet(std::istream& s);

    explicit
    IQIndexSet(const IQIndexSet& other);

    //
    // Accessor Methods
    //

    int
    r() const { return index_.size(); }

    const IQIndex&
    index(int j) const { return GET(index_,j-1); }

    int
    m(int j) const { return GET(index_,j-1).m(); }

    typedef std::vector<IQIndex>::const_iterator 
    index_it;

    //Can be used for iteration over Indices in a Foreach loop
    //e.g. Foreach(const IQIndex& I, t.index() ) { ... }
    const std::pair<index_it,index_it> 
    index() const  
        { return std::make_pair(index_.begin(),index_.end()); }

    index_it
    begin() const { return index_.begin(); }

    index_it
    end() const { return index_.end(); }


    Real
    uniqueReal() const { return ur_; }

    //
    // IQIndex Analysis
    //

    const IQIndex&
    findtype(IndexType type) const;

    bool 
    hastype(IndexType type) const;

    int 
    findindex(const IQIndex& I) const;

    const IQIndex&
    finddir(Arrow dir) const;

    bool 
    hasCommonIndex(const IQIndexSet& other) const;
    
    bool 
    hasindex(const IQIndex& I) const;

    bool 
    hasType(IndexType t) const;

    int
    minM() const;

    int
    maxM() const;

    //
    // Primelevel Methods
    //

    void 
    prime(IndexType type, int inc = 1);

    void 
    prime(const IQIndex& I, int inc = 1);

    void 
    prime(const IQIndex& I, const IQIndex& J);

    void 
    noprime(IndexType type = All);

    void 
    noprime(const IQIndex& I);

    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    void 
    mapprimeind(const IQIndex& I, int plevold, int plevnew, 
                IndexType type = All);


    // Increments primelevel of ALL copies of
    // an IQIndex I, regardless of original primelevel.
    // For example, if the original indices are 
    // I, I', J, calling indIncPrime(I,2)
    // results in I'', I''', J.
    void
    indIncAllPrime(const IQIndex& I, int inc);

    //
    // Methods for Manipulating IQIndexSets
    //
    // Warning: these can overwrite other
    // Indices if not used properly
    //

    void 
    mapindex(const IQIndex& i1, const IQIndex& i2);

    void 
    addindex(const IQIndex& I);

    void 
    removeindex(int j);

    void
    setUniqueReal();

    void
    swap(IQIndexSet& other);

    void
    swapInds(std::vector<IQIndex>& newinds);

    void
    clear();

    //
    // Other Methods
    //

    void
    conj();

    void
    conj(const IQIndex& I);

    // Checks if the IQIndexSet is properly
    // formed (e.g. no IQIndex appears more than once),
    // if not, throws an ITError.
    void
    check() const;

    void
    write(std::ostream& s) const;

    friend std::ostream&
    operator<<(std::ostream& s, const IQIndexSet& is);

    friend void 
    intrusive_ptr_add_ref(IQIndexSet* p);

    friend void 
    intrusive_ptr_release(IQIndexSet* p);

    int 
    count() const { return numref; }

    static IQIndexSet* Null()
        {
        //Set initial numref to 1000
        static IQIndexSet Null_(1000);
#ifdef DEBUG
        if(Null_.numref < 500)
            Error("Null_.numref too low");
#endif
        return &Null_;
        }

    private:

    //////////
    //
    // Data Members
    //

    std::vector<IQIndex> index_;

    Real ur_;

    mutable unsigned int 
    numref;

    //
    /////////

    explicit
    IQIndexSet(int init_numref);

    void 
    operator=(const IQIndexSet&);
    ~IQIndexSet() { }

    };


#endif
