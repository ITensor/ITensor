//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTENSOR_H
#define __ITENSOR_IQTENSOR_H
#include "iqindex.h"
#include "itensor.h"

namespace itensor {

ITensor 
toITensor(const IQTensor& T) const;

//
// IQTensor
//

class IQTensor
    {
    public:

    using IndexT = IQIndex;
    using IndexValT = IQIndexVal;
    using storage_ptr = PData;

    //Constructors --------------------------------------------------

    IQTensor();

    //Construct rank 0 IQTensor (scalar), value set to val
    //If val.imag()==0, only Real storage will be used
    explicit 
    IQTensor(Complex val);

    //Construct rank 1 IQTensor, all elements set to zero
    explicit 
    IQTensor(const IQIndex& i1);

    //Construct rank 2 IQTensor, all elements set to zero
    IQTensor(const IQIndex& i1,
             const IQIndex& i2);

    //Construct rank n IQTensor, all elements set to zero
    template<typename... IQIndices>
    IQTensor(const IQIndex& i1,
             const IQIndex& i2,
             const IQIndices&... rest);

    //Construct IQTensor with IQIndices given by vector iqinds
    explicit 
    IQTensor(std::vector<IQIndex>&& iqinds);

    //
    // IQIndexVal IQTensor Constructors
    //
    // Given a set of IQIndexVals
    // iv1 = (I1,n1), iv2 = (I2,n2), iv3 = (I3,n3), ...
    // construct an IQTensor T such that
    // T(I1(n1),I2(n2),I3(n3),...) == 1
    //
    template <typename... IQIVals>
    explicit
    IQTensor(const IQIndexVal& iv1,
             const IQIVals&... rest);

    //Accessor Methods ------------------------------------------

    //Rank of this IQTensor (number of IQIndices)
    int 
    r() const;

    //true if IQTensor has no blocks
    bool 
    empty() const;

    //IQTensor evaluates to false if it is default constructed
    explicit operator bool() const { return valid(); }

    //false if IQTensor is default constructed
    bool 
    valid() const;

    bool
    isComplex() const;

    //Returns object containing ITensor blocks
    //The ITensors can be iterated over using a Foreach
    //For example, given an IQTensor T,
    //Foreach(const ITensor& t, T.blocks()) { ... }
    const IQTDat&
    blocks() const { return *d_; }
    
    const IndexSet<IQIndex>& 
    inds() const { return is_; }


    //----------------------------------------------------
    //IQTensor operators

    //
    // Contracting product
    //
    IQTensor 
    operator*(IQTensor other) const 
        { other *= *this; return other; }

    IQTensor& 
    operator*=(const IQTensor& other);

    //
    // Addition and subtraction
    //
    IQTensor& 
    operator+=(const IQTensor& o);

    IQTensor 
    operator+(const IQTensor& o) const 
        { IQTensor res(*this); res += o; return res; }

    IQTensor&
    operator-=(const IQTensor& o)
        { 
        if(this == &o) { operator*=(0); return *this; }
        IQTensor oth(o);
        oth *= -1;
        return operator+=(oth);
        }

    IQTensor 
    operator-(const IQTensor& o) const 
        { IQTensor res(*this); res -= o; return res; }

    //
    // Multiplication by a scalar
    //
    IQTensor& 
    operator*=(Real fac);

    IQTensor& 
    operator/=(Real fac);

    IQTensor
    operator-() const { IQTensor T(*this); T *= -1; return T; }

    IQTensor& 
    operator*=(const LogNumber& lgnum);

    IQTensor& 
    operator*=(Complex z);
    
    //
    // Multiplication by an IQIndexVal
    //
    IQTensor& 
    operator*=(const IQIndexVal& iv) { return operator*=(IQTensor(iv)); }

    //Automatic conversion to ITensor
    operator 
    ITensor() const { return toITensor(); }

    //Inserts an ITensor block or adds it to
    //existing one if already present and QNs match
    IQTensor& 
    operator+=(const ITensor& block);

    //Like operator+=(ITensor) but
    //demands that the block is zero/absent
    //before inserting
    void insert(const ITensor& block);

    //----------------------------------------------------
    //IQTensor: element access

    template <typename... IQIndexVals>
    Real
    real(IQIndexVals&&... ivs) const;

    template <typename... IQIndexVals>
    Complex
    cplx(IQIndexVals&&... ivs) const;

    //Set element at location given by collection
    //of IQIndexVals. Will not switch storage
    //from Real to Complex unless val.imag()!=0 
    template<typename... IQIndexVals>
    void
    set(Complex val, const IQIndexVals&... ivs);

    //----------------------------------------------------
    //IQTensor: prime methods

    IQTensor& 
    noprime(IndexType type = All);

    IQTensor& 
    noprime(const IQIndex& I);

    IQTensor& 
    prime(int inc = 1) { prime(All,inc); return *this; }

    IQTensor& 
    prime(IndexType type, int inc = 1);

    IQTensor& 
    prime(const IQIndex& I, int inc = 1);

    //no need to keep prime level small
    IQTensor& 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //----------------------------------------------------
    //IQTensor miscellaneous methods

    IQTensor&
    takeRealPart();

    IQTensor&
    takeImagPart();

    void 
    scaleTo(const LogNumber& newscale);
    void 
    scaleTo(Real newscale) { scaleTo(LogNumber(newscale)); }

    //Take complex conjugate; do not reverse IQIndex arrows
    IQTensor& 
    conj();

    //Take complex conjugate and reverse IQIndex arrows
    IQTensor& 
    dag();

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    //
    // Deprecated methods
    //

    const IndexSet<IQIndex>& 
    indices() const { return inds(); }

    private:

    /////////////
    IQIndexSet is_;
    storage_ptr store_;
    /////////////////

    void 
    solo();

    void 
    scaleOutNorm();

    }; //class IQTensor


IQTensor inline
operator*(IQTensor T, Real fac) {  T *= fac; return T; }

IQTensor inline
operator*(Real fac, IQTensor T) { T *= fac; return T; }

IQTensor inline
operator/(IQTensor T, Real fac) {  T /= fac; return T; }

IQTensor inline
operator/(Real fac, IQTensor T) { T /= fac; return T; }

IQTensor inline
operator*(IQTensor T, const LogNumber& fac) {  T *= fac; return T; }

IQTensor inline
operator*(const LogNumber& fac, IQTensor T) { T *= fac; return T; }

IQTensor inline
operator*(IQTensor T, Complex fac) {  T *= fac; return T; }

IQTensor inline
operator*(Complex fac, IQTensor T) { T *= fac; return T; }

IQTensor inline
operator*(IQTensor T, const IQIndexVal& iv) { T *= iv; return T; }

IQTensor inline
operator*(const IQIndexVal& iv, const IQTensor& T) { return IQTensor(iv) * T; }


//
// Multiplication by an IndexVal
// Result is an ITensor
//
ITensor inline
operator*(const IQTensor& T, const IndexVal& iv)
    { 
    ITensor res = T.toITensor(); 
    return res *= iv; 
    }

ITensor inline
operator*(const IndexVal& iv, const IQTensor& T) 
    { 
    return ITensor(iv) * T.toITensor(); 
    }

//Take complex conjugate of IQTensor res,
//but do not reverse IQIndex arrows
IQTensor inline
conj(IQTensor res) { res.conj(); return res; }

//Compute complex conjugate of IQTensor res,
//and reverse IQIndex arrows
IQTensor inline
dag(IQTensor res) { res.dag(); return res; }

//Compute divergence of IQTensor T
//
//If DEBUG defined and all blocks do not have
//the same divergence, throws an exception
//(since IQTensor is not correctly constructed).
QN 
div(const IQTensor& T, const Args& args = Global::args());

const IQIndex&
findIQInd(const IQTensor& T, const Index& i);

QN
qn(const IQTensor& T, const Index& i);

Arrow
dir(const IQTensor& T, const Index& i);

Arrow
dir(const IQTensor& T, const IQIndex& i);

//Returns true if T is exactly zero.
//
//If passed the argument Args("Fast",true),
//only performs fast operations such as checking
//whether T contains any blocks, but skips computing
//the norm of the blocks.
//This can cause the return value to be true even
//if T is actually zero.
bool
isZero(const IQTensor& T, const Args& args = Global::args());

std::ostream& 
operator<<(std::ostream & s, const IQTensor &t);

}; //namespace itensor

#endif
