//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTENSOR_H
#define __ITENSOR_IQTENSOR_H
#include "iqindex.h"
#include "itensor.h"
#include "iqtdata_functions.h"

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
    valid() const { return bool(store_); }

    const IQIndexSet& 
    inds() const { return is_; }

    // Contracting product
    IQTensor& 
    operator*=(const IQTensor& other);

    // Addition and subtraction
    IQTensor& 
    operator+=(const IQTensor& o);

    IQTensor&
    operator-=(const IQTensor& o);

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

    //Add ITensor 'block' to the matching block
    //of this IQTensor. The quantum-number flux
    //of the ITensor block must be compatible
    //with this IQTensor.
    IQTensor& 
    addBlock(const ITensor& block);

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

    //Set all elements to z. If z.imag()==0
    //(such as if z is automatically converted from a Real)
    //then storage will be real only.
    IQTensor&
    fill(Complex z);

    //Call a function of the form f()->val once
    //for each element, assign result to each element.
    template <typename Func>
    IQTensor&
    generate(Func&& f);

    //Apply a function of the form f(x)->y
    //to each element x, replacing it with y
    template <typename Func>
    IQTensor&
    apply(Func&& f);

    //Apply a function of the form f(x)->void
    //to each element x.
    template <typename Func>
    const IQTensor&
    visit(Func&& f) const;

    IQTensor&
    takeRealPart();

    IQTensor&
    takeImagPart();

    const ITData&
    data() const { return *store_; }

    const LogNumber&
    scale() const { return scale_; }
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
    storage_ptr store_;
    QN div_;
    LogNumber scale_;
    /////////////////

    void 
    scaleOutNorm();

    friend QN div(const IQTensor&);

    }; //class IQTensor

template<typename... IQIndices>
IQTensor::
IQTensor(const QN& q,
         const IQIndex& i1,
         const IQIndices&... rest)
    : 
    is_(i1,rest...),
    store_(make_shared<IQTData<Real>>(is_,q)),
    div_(q),
    scale_(1.)
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
    for(const auto& iv : ivs) div_ += iv.qn()*iv.index.dir();
    store_ = make_shared<IQTData<Real>>(is_,div_);
    set(1.,iv1,rest...);
    }


template <typename... IQIndexVals>
Complex IQTensor::
cplx(IQIndexVals&&... ivs) const
    {
    constexpr auto size = sizeof...(ivs);
#ifdef DEBUG
    if(!*this) Error("IQTensor is default constructed");
    if(size > r()) Error("Too many IQIndexVals passed to real/cplx");
#endif
    std::array<IQIndexVal,size> vals{{static_cast<IQIndexVal>(ivs)...}};
    std::array<long,size> inds;
    detail::permute_map(is_,vals,inds,[](const IQIndexVal& iv) { return iv.i-1; });
    Complex z = applyFunc<IQGetElt<Complex,size>>(store_,{is_,inds});
	try {
	    return z*scale_.real(); 
	    }
	catch(const TooBigForReal& e)
	    {
	    println("too big for real in cplx(...), scale = ",scale());
	    throw e;
	    }
	catch(TooSmallForReal)
	    {
        println("warning: too small for real in cplx(...)");
	    return Complex(0.,0.);
	    }
    return Complex(NAN,NAN);
    }

template <typename... IQIndexVals>
Real IQTensor::
real(IQIndexVals&&... ivs) const
    {
    auto z = cplx(std::forward<IQIndexVals>(ivs)...);
    if(fabs(z.imag()) != 0)
        {
        printfln("element = (%.5E,%.5E)",z.real(),z.imag());
        Error("IQTensor is Complex-valued, use .cplx(...) method");
        }
    return z.real();
    }

template<typename... IQIndexVals>
void IQTensor::
set(Complex val, const IQIndexVals&... ivs)
    {
    static constexpr auto size = sizeof...(ivs);
    scaleTo(1.);
    const std::array<IQIndexVal,size> vals{{ static_cast<IQIndexVal>(ivs)...}};
    std::array<long,size> inds;
    detail::permute_map(is_,vals,inds,[](const IQIndexVal& iv) { return iv.i-1; });
    if(val.imag() == 0)
        applyFunc<IQSetEltReal<size>>(store_,{val.real(),is_,inds});
    else
        {
        Error("Set Complex element not implemented");
        //applyFunc<IQSetEltComplex<size>>(store_,{val,is_,bis});
        }
    }

template <typename Func>
IQTensor& IQTensor::
generate(Func&& f)
    {
    scaleTo(1);
    applyFunc<GenerateIQT<decltype(f)>>(store_,{std::forward<Func>(f)});
    return *this;
    }

template <typename Func>
IQTensor& IQTensor::
apply(Func&& f)
    {
    scaleTo(1);
    applyFunc<ApplyIQT<decltype(f)>>(store_,{std::forward<Func>(f)});
    return *this;
    }

template <typename Func>
const IQTensor& IQTensor::
visit(Func&& f) const
    {
    applyFunc<VisitIQT<decltype(f)>>(store_,{std::forward<Func>(f),scale_.real0()});
    return *this;
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

bool
isComplex(const IQTensor& T);

//Take complex conjugate of IQTensor res,
//but do not reverse IQIndex arrows
IQTensor inline
conj(IQTensor res) { res.conj(); return res; }

//Compute complex conjugate of IQTensor res,
//and reverse IQIndex arrows
IQTensor inline
dag(IQTensor res) { res.dag(); return res; }

//Compute divergence of IQTensor T
QN inline
div(const IQTensor& T) 
    { 
    if(!T) Error("div(IQTensor) not defined for null IQTensor");
    return T.div_; 
    }

IQIndex
findIQInd(const IQTensor& T, const Index& i);

QN inline
qn(const IQTensor& T, const Index& i) { return qn(findIQInd(T,i),i); }

Arrow inline
dir(const IQTensor& T, const Index& i) { return findIQInd(T,i).dir(); }

Arrow
dir(const IQTensor& T, const IQIndex& i);

//Compute the norm of an IQTensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt((T*T).real()) 
//[and similar for complex case] but computed much more efficiently.
Real
norm(const IQTensor& T);

IQTensor
randomize(IQTensor T, const Args& args = Global::args());

template <typename... Params>
IQTensor
randIQT(Params&&... params)
    {
    return randomize(IQTensor(std::forward<Params>(params)...));
    }


std::ostream& 
operator<<(std::ostream & s, const IQTensor &t);

}; //namespace itensor

#endif
