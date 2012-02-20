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
    findtype(IndexType t) const;

    bool 
    findtype(IndexType t, IQIndex& I) const;

    int 
    findindex(const IQIndex& I) const;

    bool 
    has_common_index(const IQTSparse& other) const;
    
    bool 
    hasindex(const IQIndex& I) const;

    bool 
    notin(const IQIndex& I) const { return !hasindex(I); }

    //
    // Primelevel Methods 
    //

    void 
    noprime(PrimeType p = primeBoth);

    void 
    doprime(PrimeType pt, int inc = 1);

    void 
    primeall() { doprime(primeBoth,1); }

    void 
    primesite() { doprime(primeSite,1); }

    void 
    primelink() { doprime(primeLink,1); }

    void 
    mapprime(int plevold, int plevnew, PrimeType pt = primeBoth);

    void 
    mapprimeind(const IQIndex& I, int plevold, int plevnew, 
                PrimeType pt = primeBoth);

    void 
    primeind(const IQIndex& I, int inc = 1);

    void 
    noprimeind(const IQIndex& I);

    //
    // Other Methods
    //

    template <typename Callable> void
    mapElems(const Callable& f); 

    void
    conj();

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
