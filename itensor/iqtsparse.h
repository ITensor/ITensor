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

    IQTSparse(const IQIndex& i1);

    IQTSparse(const IQIndex& i1, const IQIndex& i2);

    IQTSparse(const IQIndex& i1, const IQIndex& i2, 
                  const IQIndex& i3);

    IQTSparse(const IQIndex& i1, const IQIndex& i2, 
                  const IQIndex& i3, const IQIndex& i4);

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
    isNull() const { return (d_ == 0 || is_ == 0); }

    bool
    isNotNull() const { return (d_ != 0 && is_ != 0); }

    bool 
    isComplex() const { return hasindex(IQIndex::IndReIm()); }

    bool 
    isNotComplex() const { return !hasindex(IQIndex::IndReIm()); }

    bool
    isDiag() const { return true; }

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

    //
    // IQIndex Methods
    //

    IQIndex 
    findtype(IndexType t) const { return is_->findtype(t); }

    bool 
    findtype(IndexType t, IQIndex& I) const { return is_->findtype(t,I); }

    int 
    findindex(const IQIndex& I) const { return is_->findindex(I); }

    bool 
    has_common_index(const IQTSparse& other) const
        { return is_->has_common_index(*other.is_); }
    
    bool 
    hasindex(const IQIndex& I) const { return is_->hasindex(I); }

    bool 
    notin(const IQIndex& I) const { return !hasindex(I); }

    void 
    mapindex(const IQIndex& i1, const IQIndex& i2) { is_->mapindex(i1,i2); }

    //
    // Primelevel Methods 
    //

    void 
    noprime(PrimeType p = primeBoth) { is_->noprime(p); }

    void 
    doprime(PrimeType pt, int inc = 1) { is_->doprime(pt,inc); }

    void 
    primeall() { doprime(primeBoth,1); }

    void 
    primesite() { doprime(primeSite,1); }

    void 
    primelink() { doprime(primeLink,1); }

    void 
    mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
        { is_->mapprime(plevold,plevnew,pt); }

    void 
    mapprimeind(const IQIndex& I, int plevold, int plevnew, 
                PrimeType pt = primeBoth)
        { is_->mapprimeind(I,plevold,plevnew,pt); }

    void 
    primeind(const IQIndex& I, int inc = 1)
        { mapindex(I,primed(I,inc)); }

    void 
    primeind(const IQIndex& I, const IQIndex& J) { is_->primeind(I,J); }

    void 
    noprimeind(const IQIndex& I) { is_->noprimeind(I); }

    //
    // Other Methods
    //

    Real
    norm() const;

    void 
    scaleOutNorm() const;

    void 
    scaleTo(const LogNumber& newscale) const;

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    friend std::ostream&
    operator<<(std::ostream & s, const IQTSparse & t);

    typedef IQIndex 
    IQIndexT;

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

    IQTSDat(std::istream& s);

    //
    // Accessors
    //

    const_iterator
    begin() const { return its_.begin(); }
    iterator
    begin() { return its_.begin(); }

    const_iterator
    end() const { return its_.end(); }
    iterator
    end() { return its_.end(); }

    void
    insert_add(const ITSparse& s);

    void
    clear();

    //
    // Other Methods
    //

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    ENABLE_INTRUSIVE_PTR(IQTSDat)

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

    void
    uninit_rmap() const;

    void
    init_rmap() const;

    };

void 
product(const IQTSparse& S, const IQTensor& T, IQTensor& res);

#endif
