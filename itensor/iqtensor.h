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
// IQTensor related functions
//

//Specialization of ITensorT::dag()
template<>
IQTensor& IQTensor::dag();

// Contract with IndexVal
// If iv = (J,n), Index J is fixed to it's nth
// value and rank decreases by 1
// (similar to summing against a Kronecker
// delta tensor \delta_{J,n})
inline IQTensor& 
operator*=(IQTensor& T, const IQIndexVal& iv) { return T *= IQTensor(iv); } 
IQTensor inline
operator*(IQTensor T, const IQIndexVal& iv) { T *= iv; return T; }
IQTensor inline
operator*(const IQIndexVal& iv, const IQTensor& T) { return IQTensor(iv) * T; }


//Add ITensor to corresponding block of IQTensor
IQTensor& 
operator+=(IQTensor& A, const ITensor& B);

ITensor 
toITensor(const IQTensor& T);

template<> inline
IQTensor::
operator ITensor() const { return toITensor(*this); }

ITensor inline
operator*(const IQTensor& T, const ITensor& t) 
    { 
    auto TT = toITensor(T);
    TT *= t; 
    return TT; 
    }
ITensor inline
operator*(const ITensor& t, const IQTensor& T) 
    { 
    return operator*(T,t);
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
