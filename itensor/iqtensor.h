//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTENSOR_H
#define __ITENSOR_IQTENSOR_H
#include "itensor/itensor.h"
#include "itensor/iqindex.h"

namespace itensor {

//
// IQTensor
//

using IQTensor = ITensorT<IQIndex>;


template<>
IQTensor& IQTensor::conj();

template<>
IQTensor& IQTensor::dag();

template<>
void 
IQTensor::scaleTo(const LogNumber&);

//Contracting product
//All matching Index pairs automatically contracted
//Cji = \sum_{k,l} Akjl * Blki
IQTensor& 
operator*=(IQTensor& A, const IQTensor& B);

// Contract with IndexVal
// If iv = (J,n), Index J is fixed to it's nth
// value and rank decreases by 1
// (similar to summing against a Kronecker
// delta tensor \delta_{J,n})
inline IQTensor& 
operator*=(IQTensor& T, const IQIndexVal& iv) { return operator*=(T,IQTensor(iv)); } 

//Tensor addition and subtraction
//Summands must have same Indices, in any order
//Cijk = Aijk + Bkij
IQTensor& 
operator+=(IQTensor& A, const IQTensor& B);
IQTensor& 
operator-=(IQTensor& A, const IQTensor& B);

//Add ITensor to corresponding block of IQTensor
IQTensor& 
operator+=(IQTensor& A, const ITensor& B);

//TODO
//Multiplication and division by complex scalar
//IQTensor& 
//operator*=(IQTensor& T, Cplx z);
//inline IQTensor& 
//operator/=(IQTensor& T, Cplx z) { return operator*=(T,1./z); }

IQTensor inline
operator*(IQTensor A, const IQTensor& B) { A *= B; return A; }
IQTensor inline
operator*(IQTensor T, Real fac) {  T *= fac; return T; }
IQTensor inline
operator*(Real fac, IQTensor T) { T *= fac; return T; }
IQTensor inline
operator/(IQTensor T, Real fac) {  T /= fac; return T; }

//TODO
//IQTensor inline
//operator*(IQTensor T, Complex fac) {  T *= fac; return T; }
//IQTensor inline
//operator*(Complex fac, IQTensor T) { T *= fac; return T; }

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

template<> inline
IQTensor::
operator ITensor() const { return toITensor(*this); }

// Read IQTensor from binary input stream.
void inline
read(std::istream& s, IQTensor& t)
    {
    Error("IQTensor read not yet implemented");
    }

// Write IQTensor to binary output stream.
void inline
write(std::ostream& s, const IQTensor& t)
    {
    Error("IQTensor read not yet implemented");
    }

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
QN
div(const IQTensor& T);

IQTensor
combiner(std::vector<IQIndex> inds, const Args& args = Global::args());

template<typename... Inds>
IQTensor
combiner(const IQIndex& i1, const Inds&... inds)
    {
    return combiner(std::vector<IQIndex>{i1,inds...});
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

//template <typename... Inds>
//IQTensor
//randomTensor(const IQIndex& i1, Inds&&... inds)
//    {
//    return randomize(IQTensor(i1,std::forward<Inds>(inds)...));
//    }
template <typename... IndVals>
IQTensor
randomTensor(const IQIndexVal& iv1, IndVals&&... ivs)
    {
    return randomize(IQTensor(iv1,std::forward<IndVals>(ivs)...));
    }
template <typename... Inds>
IQTensor
randomTensor(const QN& q, const IQIndex& i1, Inds&&... inds)
    {
    auto is = IQIndexSet{i1,std::forward<Inds>(inds)...};
    auto dat = IQTData{is,q};
    return randomize(IQTensor(std::move(is),std::move(dat)));
    }

std::ostream& 
operator<<(std::ostream & s, const IQTensor &t);

} //namespace itensor

//See file iqtensor.ih for template/inline method implementations
#include "itensor/iqtensor.ih"

#endif
