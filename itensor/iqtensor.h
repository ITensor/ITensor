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
IQTensor& 
operator*=(IQTensor& T, IQIndexVal const& iv);
IQTensor
operator*(IQTensor T, IQIndexVal const& iv);
IQTensor
operator*(IQIndexVal const& iv, IQTensor const& T);
ITensor
operator*(IndexVal const& iv, IQTensor const& T);


// Contract ITensor with IQIndexVal
ITensor& 
operator*=(ITensor& T, IQIndexVal const& iv);
ITensor
operator*(ITensor T, IQIndexVal const& iv);
ITensor
operator*(IQIndexVal const& iv, ITensor const& T);

//Add ITensor to corresponding block of IQTensor
IQTensor& 
operator+=(IQTensor& A, ITensor const& B);

ITensor 
toITensor(IQTensor const& T);

template<> inline
IQTensor::
operator ITensor() const { return toITensor(*this); }

ITensor 
operator*(IQTensor const& T, ITensor const& t);

ITensor
operator*(ITensor const& t, IQTensor const& T);


//
// Multiplication by an IndexVal
// Result is an ITensor
//
ITensor
operator*(IQTensor const& T, IndexVal const& iv);
ITensor
operator*(IndexVal const& iv, IQTensor const& T);

//Compute flux (a.k.a. divergence) of IQTensor T
QN
flux(IQTensor const& T);
QN
div(IQTensor const& T);

IQTensor
combiner(std::vector<IQIndex> inds, Args const& args = Global::args());

template<typename... Inds>
IQTensor
combiner(IQIndex const& i1, 
         Inds const&... inds)
    {
    return combiner(std::vector<IQIndex>{i1,inds...});
    }

IQTensor
delta(IQIndex const& i1, IQIndex const& i2);

IQIndex
findIQInd(const IQTensor& T, const Index& i);

QN inline
qn(const IQTensor& T, const Index& i) { return qn(findIQInd(T,i),i); }

Arrow inline
dir(const IQTensor& T, const Index& i) { return findIQInd(T,i).dir(); }

Arrow
dir(const IQTensor& T, const IQIndex& i);

template <typename... Inds>
IQTensor
randomTensor(IQIndex const& i1, Inds&&... inds)
    {
    static_assert(stdx::false_regardless_of<Inds...>::value,"Must provide a QN or IQIndexVals to IQIndex version of randomTensor");
    return IQTensor{};
    }

template <typename... IQIndVals>
IQTensor
randomTensor(IQIndexVal const& iv1, 
             IQIndVals&&... ivs)
    {
    auto T = setElt(iv1,std::forward<IQIndVals>(ivs)...);
    try {
        return random(T);
        }
    catch(ITError const& e)
        {
        Error("Cannot randomize IQTensor, possibly incompatible flux");
        }
    return T;
    }
template <typename... Inds>
IQTensor
randomTensor(QN const& q, IQIndex const& i1, Inds &&... inds)
    {
    auto is = IQIndexSet{i1,std::forward<Inds>(inds)...};
    auto dat = QDenseReal{is,q};
    auto T = IQTensor(std::move(is),std::move(dat));
    try {
        return random(T);
        }
    catch(ITError const& e)
        {
        Error("Cannot randomize IQTensor, possibly incompatible flux");
        }
    return T;
    }

template <typename... VArgs>
IQTensor
randomTensorC(QN const& q, VArgs&&... vargs)
    {
    auto T = randomTensor(q,std::forward<VArgs>(vargs)...);
    return T+1_i*random(T);
    }

//mixedIQTensor constructs
//an IQTensor with MixedQN storage
//allowing setting elements in
//multiple QN sectors.
//This is useful if creating an IQTensor
//whose only purpose is to be converted
//to an ITensor.
template<typename... Inds>
IQTensor
mixedIQTensor(IQIndex const& i1, 
              Inds const&... inds);

std::ostream& 
operator<<(std::ostream & s, IQTensor const& t);

} //namespace itensor

//See file iqtensor.ih for template/inline method implementations
#include "itensor/iqtensor.ih"

#endif
