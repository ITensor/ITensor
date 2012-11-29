//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTSPARSE_H
#define __ITENSOR_IQTSPARSE_H
#include "itsparse.h"
#include "iqtensor.h"

class IQTSDat;

//
// IQTSparse
//
class IQTSparse
    {
    public:

    //
    // Constructors
    //

    IQTSparse();

    explicit
    IQTSparse(const IQIndex& i1);

    IQTSparse(const IQIndex& i1, const IQIndex& i2);

    //Set diagonal to a constant d. If d == 1 makes a Kronecker delta
    IQTSparse(const IQIndex& i1, const IQIndex& i2, Real d);

    IQTSparse(const IQIndex& i1, const IQIndex& i2, 
                  const IQIndex& i3);

    IQTSparse(const IQIndex& i1, const IQIndex& i2, 
                  const IQIndex& i3, const IQIndex& i4);

    explicit
    IQTSparse(std::istream& s) { read(s); }

    //
    // Accessor Methods
    //

    const IQIndex& 
    index(int j) const { return is_->index(j); }

    int 
    r() const { return is_->r(); }

    int 
    m(int j) const { return is_->m(j); }

    //uniqueReal depends on indices only, unordered:
    Real 
    uniqueReal() const { return is_->uniqueReal(); } 

    bool
    isNull() const;

    bool
    isNotNull() const;

    bool 
    isComplex() const { return hasindex(IQIndex::IndReIm()); }

    bool 
    isNotComplex() const { return !hasindex(IQIndex::IndReIm()); }

    bool
    isDiag() const { return true; }

    const IQTSDat&
    blocks() const { return *d_; }

    //
    // Operators
    //

    // Contracting product with IQTensors

    IQTensor
    operator*(const IQTensor& T) const
        { IQTensor res; product(*this,T,res); return res; }

    friend inline IQTensor&
    operator*=(IQTensor& T, const IQTSparse& S)
        { IQTensor res; product(S,T,res); T.swap(res); return T; }

    IQTensor friend inline
    operator*(const IQTensor& T, const IQTSparse& S)
        { IQTensor res; product(S,T,res); return res; }

    // Addition with ITSparse

    IQTSparse&
    operator+=(const ITSparse& s);

    // Addition with IQTSparse

    IQTSparse&
    operator+=(const IQTSparse& s);

    // Multiplication and division by a scalar

    IQTSparse& 
    operator*=(Real fac);

    IQTSparse 
    operator*(Real fac) const 
        { IQTSparse res(*this); res *= fac; return res; }

    friend inline IQTSparse 
    operator*(Real fac, IQTSparse s) 
        { return (s *= fac); }

    IQTSparse& 
    operator/=(Real fac) { operator*=(1./fac); return *this; }

    IQTSparse 
    operator/(Real fac) const 
        { IQTSparse res(*this); res /= fac; return res; }

    friend inline IQTSparse 
    operator/(Real fac, IQTSparse s) { return (s /= fac); }

    IQTSparse& 
    operator*=(const LogNumber& fac);

    IQTSparse 
    operator*(const LogNumber& fac) const 
        { IQTSparse res(*this); res *= fac; return res; }

    friend inline IQTSparse 
    operator*(const LogNumber& fac, IQTSparse s) 
        { return (s *= fac); }

    //Real& 
    //operator()(const IQIndexVal& iv1, 
    //           const IQIndexVal& iv2 = IQIndexVal::Null(), 
	//           const IQIndexVal& iv3 = IQIndexVal::Null(), 
    //           const IQIndexVal& iv4 = IQIndexVal::Null(), 
	//           const IQIndexVal& iv5 = IQIndexVal::Null(), 
    //           const IQIndexVal& iv6 = IQIndexVal::Null(),
	//           const IQIndexVal& iv7 = IQIndexVal::Null(), 
    //           const IQIndexVal& iv8 = IQIndexVal::Null());

    //
    // IQIndex Methods
    //

    IQIndex 
    findtype(IndexType t) const;

    int 
    findindex(const IQIndex& I) const;

    bool 
    has_common_index(const IQTSparse& other) const;
    
    bool 
    hasindex(const IQIndex& I) const;

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
    mapprimeind(const IQIndex& I, int plevold, int plevnew, 
                IndexType type = All);

    void 
    primeind(const IQIndex& I, int inc = 1);

    void 
    noprimeind(const IQIndex& I);

    friend inline IQTSparse
    primed(IQTSparse S, int inc = 1)
        { S.doprime(All,inc); return S; }

    friend inline IQTSparse
    primesite(IQTSparse S, int inc = 1)
        { S.doprime(Site,inc); return S; }

    friend inline IQTSparse
    primelink(IQTSparse S, int inc = 1)
        { S.doprime(Link,inc); return S; }

    friend inline IQTSparse
    primeind(IQTSparse S, const IQIndex& I, int inc = 1)
        { S.primeind(I,inc); return S; }

    friend inline IQTSparse
    deprimed(IQTSparse S)
        { S.noprime(); return S; }

    //
    // Other Methods
    //

    template <typename Callable> void
    mapElems(const Callable& f); 

    void
    conj();

    void
    pseudoInvert(Real cutoff = 0);

    Real
    norm() const;

    void 
    scaleOutNorm() const;

    void 
    scaleTo(const LogNumber& newscale) const;

    void 
    print(std::string name = "",Printdat pdat = HideData) const;

    void 
    printIndices(const std::string& name = "") const;

    inline void 
    printIndices(const boost::format& fname) const
        { printIndices(fname.str()); }


    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    friend std::ostream&
    operator<<(std::ostream & s, const IQTSparse & t);

    typedef IQIndex 
    IndexT;

    typedef IQIndexVal 
    IQIndexValT;

    typedef Combiner 
    CombinerT;

    private:

    //////////////
    //
    // Data members
    //

    boost::intrusive_ptr<IQIndexSet> is_;

    boost::intrusive_ptr<IQTSDat> d_;

    //
    //////////////

    void
    soloDat();

    void
    soloIndex();

    void
    solo();


    IQTSDat&
    ncblocks() { return *d_; }

    friend void 
    product(const IQTSparse& S, const IQTensor& T, IQTensor& res);

    }; // class IQTSparse

//
// IQTSDat (Storage for IQTSparse)
//

class IQTSDat
    {
    public:

    typedef std::list<ITSparse>
    StorageT;

    typedef StorageT::const_iterator
    const_iterator;

    typedef StorageT::iterator
    iterator;

    //
    // Constructors
    //

    IQTSDat();

    explicit
    IQTSDat(const IQTSDat& other);

    explicit
    IQTSDat(std::istream& s);

    //
    // Accessors
    //

    const_iterator
    begin() const { return its_.begin(); }

    iterator
    begin() { uninit_rmap(); return its_.begin(); }

    const_iterator
    end() const { return its_.end(); }

    iterator
    end() { uninit_rmap(); return its_.end(); }

    void
    insert_add(const ITSparse& s);

    void
    clear();

    //
    // Other Methods
    //

    void 
    scaleTo(const LogNumber& newscale) const;

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    friend void 
    intrusive_ptr_add_ref(IQTSDat* p);

    friend void 
    intrusive_ptr_release(IQTSDat* p);

    int 
    count() const { return numref; }

    static IQTSDat* Null()
        {
        //Set initial numref to 1000, stack allocated
        static IQTSDat Null_(1000);
#ifdef DEBUG
        if(Null_.numref < 500)
            Error("Null_.numref too low");
#endif
        return &Null_;
        }

    private:

    //////////////
    //
    // Data members
    //

    mutable unsigned int numref;

    mutable bool init;

    mutable StorageT its_;

    mutable std::map<ApproxReal,iterator>
    rmap;

    //
    //////////////

    explicit
    IQTSDat(int init_numref);

    void
    uninit_rmap() const;

    void
    init_rmap() const;

    };

void 
product(const IQTSparse& S, const IQTensor& T, IQTensor& res);

template <typename Callable> void IQTSparse::
mapElems(const Callable& f)
    {
    soloDat();

    Foreach(ITSparse& s, *(d_))
        { 
        s.mapElems(f);
        }
    }


#endif
