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

    //Set diagonal to a vector D.
    IQTSparse(const IQIndex& i1, const IQIndex& i2, const VectorRef& D);

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
    isDiag() const { return true; }

    const IQTSDat&
    blocks() const { return *d_; }

    const IndexSet<IQIndex>& 
    indices() const { return *is_; }

    bool
    isComplex() const;

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
    // Primelevel Methods 
    //

    IQTSparse& 
    prime(int inc = 1) { prime(All,inc); return *this; }

    IQTSparse& 
    prime(const IQIndex& I, int inc = 1);

    IQTSparse& 
    prime(IndexType type, int inc = 1);

    IQTSparse& 
    noprime(IndexType type = All);

    IQTSparse& 
    noprime(const IQIndex& I);

    IQTSparse& 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //
    // Other Methods
    //

    template <typename Callable> 
    IQTSparse&
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

    boost::shared_ptr<IndexSet<IQIndex> > is_;

    boost::shared_ptr<IQTSDat> d_;

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
    makeCopyOf(const IQTSDat& other);

    void 
    scaleTo(const LogNumber& newscale) const;

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    static const boost::shared_ptr<IQTSDat>& Null();

    private:

    //////////////
    //
    // Data members
    //

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

template <typename Callable> 
IQTSparse& IQTSparse::
mapElems(const Callable& f)
    {
    soloDat();

    Foreach(ITSparse& s, *(d_))
        { 
        s.mapElems(f);
        }
    return *this;
    }


#endif
