//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTENSOR_H
#define __ITENSOR_IQTENSOR_H
#include "iqindex.h"
#include "itensor.h"
#include "itdata/iqtdense.h"

namespace itensor {

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

    IQTensor() { }

    //Construct rank 0 IQTensor (scalar), value set to val
    //If val.imag()==0, only Real storage will be used
    explicit 
    IQTensor(Complex val);

    //Construct rank n IQTensor T with divergence
    //div(T)==q, all elements set to zero
    template<typename... IQIndices>
    explicit
    IQTensor(const QN& q,
             const IQIndex& i1,
             const IQIndices&... rest);

    //Construct IQTensor with IQIndices given by vector iqinds
    explicit 
    IQTensor(const QN& q, 
             std::vector<IQIndex>&& iqinds);

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
    r() const { return is_.r(); }

    //IQTensor evaluates to false if it is default constructed
    explicit operator bool() const { return valid(); }

    //false if IQTensor is default constructed
    bool 
    valid() const;

    const IQIndexSet& 
    inds() const { return is_; }

    // Contracting product
    IQTensor& 
    operator*=(const IQTensor& other);

    // Addition and subtraction
    IQTensor& 
    operator+=(const IQTensor& o);

    IQTensor&
    operator-=(const IQTensor& o)
        { 
        if(this == &o) { operator*=(0); return *this; }
        IQTensor oth(o);
        oth *= -1;
        return operator+=(oth);
        }

    //
    // Multiplication by a scalar
    //
    IQTensor& 
    operator*=(Real fac);

    IQTensor& 
    operator/=(Real fac);

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
    operator ITensor() const;

    //Inserts an ITensor block or adds it to
    //existing one if already present and QNs match
    IQTensor& 
    operator+=(const ITensor& block);

    //Like operator+=(ITensor) but
    //demands that the block is zero/absent
    //before inserting
    void 
    insert(const ITensor& block);

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
    noprime(IndexType type = All) { is_.noprime(type); return *this; }

    IQTensor& 
    noprime(const IQIndex& I) { is_.noprime(I); return *this; }

    IQTensor& 
    prime(int inc = 1) { prime(All,inc); return *this; }

    IQTensor& 
    prime(IndexType type, int inc = 1) { is_.prime(type,inc); return *this; }

    IQTensor& 
    prime(const IQIndex& I, int inc = 1) { is_.prime(I,inc); return *this; }

    //no need to keep prime level small
    IQTensor& 
    mapprime(int plevold, int plevnew, IndexType type = All) 
        { is_.mapprime(plevold,plevnew,type); return *this; }

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

    const IQIndexSet& 
    indices() const { return inds(); }

    private:

    /////////////////
    IQIndexSet is_;
    LogNumber scale_;
    storage_ptr store_;
    /////////////////

    void 
    solo();

    void 
    scaleOutNorm();

    void
    initDense(const QN& Q);

    }; //class IQTensor

template<typename... IQIndices>
IQTensor::
IQTensor(const QN& q,
         const IQIndex& i1,
         const IQIndices&... rest)
    : 
    is_(i1,rest...),
    scale_(1.),
    store_(make_shared<IQTDense<Real>>(is_,q))
    { }

template <typename... IQIVals>
IQTensor::
IQTensor(const IQIndexVal& iv1,
         const IQIVals&... rest)
    :
    scale_(1.)
    { 
    const size_t size = 1+sizeof...(rest);
    auto ivs = std::array<IQIndexVal,size>{{iv1,rest...}};
    std::array<IQIndex,size> inds;
    for(size_t j = 0; j < size; ++j) inds[j] = ivs[j].index;
    is_ = IQIndexSet(inds);
    QN q;
    for(const auto& iv : ivs) q += iv.qn()*iv.index.dir();
    store_ = make_shared<IQTDense<Real>>(is_,q);
    set(1.,iv1,rest...);
    }


IQTensor inline
operator*(IQTensor A, const IQTensor& B) { A *= B; return A; }
IQTensor inline
operator*(IQTensor T, Real fac) {  T *= fac; return T; }
IQTensor inline
operator*(Real fac, IQTensor T) { T *= fac; return T; }
IQTensor inline
operator/(IQTensor T, Real fac) {  T /= fac; return T; }
IQTensor inline
operator*(IQTensor T, Complex fac) {  T *= fac; return T; }
IQTensor inline
operator*(Complex fac, IQTensor T) { T *= fac; return T; }
IQTensor inline
operator*(IQTensor T, const LogNumber& fac) {  T *= fac; return T; }
IQTensor inline
operator*(const LogNumber& fac, IQTensor T) { T *= fac; return T; }
IQTensor inline
operator*(IQTensor T, const IQIndexVal& iv) { T *= iv; return T; }
IQTensor inline
operator*(const IQIndexVal& iv, const IQTensor& T) { return IQTensor(iv) * T; }
IQTensor inline
operator+(IQTensor A, const IQTensor& B) { A += B; return A; }
IQTensor inline
operator-(IQTensor A, const IQTensor& B) { A -= B; return A; }
IQTensor inline
operator-(IQTensor T) { T *= -1; return T; }

ITensor 
toITensor(const IQTensor& T);

inline IQTensor::
operator ITensor() const { return toITensor(*this); }

//
// Multiplication by an IndexVal
// Result is an ITensor
//
ITensor inline
operator*(const IQTensor& T, const IndexVal& iv)
    { 
    return toITensor(T)*iv; 
    }

ITensor inline
operator*(const IndexVal& iv, const IQTensor& T) 
    { 
    return ITensor(iv) * toITensor(T); 
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

IQIndex
findIQInd(const IQTensor& T, const Index& i);

QN inline
qn(const IQTensor& T, const Index& i) { return qn(findIQInd(T,i),i); }

Arrow inline
dir(const IQTensor& T, const Index& i) { return findIQInd(T,i).dir(); }

Arrow
dir(const IQTensor& T, const IQIndex& i);

std::ostream& 
operator<<(std::ostream & s, const IQTensor &t);

}; //namespace itensor

#endif
