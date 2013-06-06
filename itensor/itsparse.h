//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITSPARSE_H
#define __ITENSOR_ITSPARSE_H
#include "itensor.h"

//
// ITSparse
//
class ITSparse
    {
    public:

    //
    // Constructors
    //

    ITSparse();

    explicit
    ITSparse(const Index& i1);

    ITSparse(const Index& i1, Real d);

    ITSparse(const Index& i1, const Vector& diag);

    ITSparse(const Index& i1, const Index& i2, Real d);

    ITSparse(const Index& i1, const Index& i2, const Vector& diag);

    ITSparse(const Index& i1, const Index& i2, 
                  const Index& i3, Real d);

    ITSparse(const Index& i1, const Index& i2, 
                  const Index& i3, const Vector& diag);


    ITSparse(const Index& i1, const Index& i2, 
                  const Index& i3, const Index& i4, Real d);

    explicit
    ITSparse(std::istream& s) { read(s); }

    //
    // Accessor Methods
    //

    const Index& 
    index(int j) const { return is_.index(j); }

    int 
    r() const { return is_.r(); }

    int 
    rn() const { return is_.rn(); }

    int 
    m(int j) const { return is_.m(j); }

    //uniqueReal depends on indices only, unordered:
    Real 
    uniqueReal() const { return is_.uniqueReal(); } 

    const LogNumber&
    scale() const { return scale_; }

    bool
    isDiag() const { return true; }

    bool
    diagAllSame() const { return diag_.Length() == 0; }

    int
    diagSize() const;

    Vector
    diag() const;

    void
    diag(VectorRef v);

    //Enables looping over Indices in a Foreach
    //e.g. Foreach(const Index& I, t.index() ) { ... }
    //const std::pair<IndexSet<Index>::index_it,IndexSet<Index>::index_it> 
    const IndexSet<Index>&
    indices() const { return is_; }

    bool
    isComplex() const { return false; }

    //
    // Operators
    //

    // Contracting product with ITensor

    ITensor
    operator*(const ITensor& T) const
        { ITensor res; product(*this,T,res); return res; }

    friend inline ITensor&
    operator*=(ITensor& T, const ITSparse& S)
        { ITensor res; product(S,T,res); T.swap(res); return T; }

    ITensor friend inline
    operator*(const ITensor& T, const ITSparse& S)
        { ITensor res; product(S,T,res); return res; }

    // Addition and subtraction

    ITSparse&
    operator+=(const ITSparse& other);

    ITSparse
    operator+(const ITSparse& other)
        { ITSparse res(*this); res += other; return res; }

    ITSparse&
    operator-=(const ITSparse& other);

    ITSparse
    operator-(const ITSparse& other)
        { ITSparse res(*this); res -= other; return res; }

    // Multiplication and division by scalars

    ITSparse& 
    operator*=(Real fac) { scale_ *= fac; return *this; }

    ITSparse 
    operator*(Real fac) const 
        { ITSparse res(*this); res *= fac; return res; }

    friend inline ITSparse 
    operator*(Real fac, ITSparse s) 
        { return (s *= fac); }

    ITSparse& 
    operator/=(Real fac) { scale_ /= fac; return *this; }

    ITSparse 
    operator/(Real fac) const 
        { ITSparse res(*this); res /= fac; return res; }

    friend inline ITSparse 
    operator/(Real fac, ITSparse s) 
        { return (s /= fac); }

    ITSparse& 
    operator*=(const LogNumber& fac) { scale_ *= fac; return *this; }

    ITSparse 
    operator*(const LogNumber& fac) const 
        { ITSparse res(*this); res *= fac; return res; }

    friend inline ITSparse 
    operator*(const LogNumber& fac, ITSparse s) 
        { return (s *= fac); }

    //
    // Index Methods
    //

    //void 
    //addindex1(const std::vector<Index>& indices) { is_.addindex1(indices); }

    //void 
    //addindex1(const Index& I) { is_.addindex1(I); }

    //
    // Primelevel Methods 
    //

    ITSparse& 
    prime(int inc = 1) { is_.prime(inc); return *this; }

    ITSparse& 
    prime(IndexType type, int inc = 1) { is_.prime(type,inc); return *this; }

    ITSparse& 
    prime(const Index& I, int inc = 1) { is_.prime(I,inc); return *this; }

    ITSparse& 
    noprime(IndexType type = All) { is_.noprime(type); return *this; }

    ITSparse& 
    noprime(const Index& I) { is_.noprime(I); return *this; }

    ITSparse& 
    mapprime(int plevold, int plevnew, IndexType type = All)
        { is_.mapprime(plevold,plevnew,type); return *this; }

    //
    // Other Methods
    //

    template <typename Callable> 
    ITSparse&
    mapElems(const Callable& f);

    void
    pseudoInvert(Real cutoff = 0);

    Real
    norm() const;

    void 
    scaleOutNorm() const;

    void 
    scaleTo(LogNumber newscale) const;

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    void
    conj() { }

    bool
    isNull() const { return (scale_ == LogNumber(0) && diag_.Length() == 0); }

    bool
    isNotNull() const { return (scale_ != LogNumber(0) || diag_.Length() != 0); }

    friend std::ostream&
    operator<<(std::ostream & s, const ITSparse & t);

    typedef Index 
    IndexT;

    typedef IndexVal 
    IndexValT;

    typedef Combiner 
    CombinerT;

    private:

    //////////////
    //
    // Data members
    //

    //diagonal elements
    mutable Vector diag_;

    IndexSet<Index> is_;

    mutable LogNumber scale_;

    //
    //////////////

    void
    _construct1(const Index& i1);

    void
    _construct2(const Index& i1, const Index& i2);

    friend void 
    product(const ITSparse& S, const ITensor& T, ITensor& res);

    }; // class ITSparse

void 
product(const ITSparse& S, const ITensor& T, ITensor& res);

template <typename Callable> 
ITSparse& ITSparse::
mapElems(const Callable& f)
    {
    scaleTo(1);
    for(int j = 1; j <= diag_.Length(); ++j)
        diag_(j) = f(diag_(j));
    return *this;
    }

#endif
