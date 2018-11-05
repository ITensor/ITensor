//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTENSOR_H
#define __ITENSOR_IQTENSOR_H
#include "itensor/itensor.h"

namespace itensor {

//Compute divergence of IQTensor T
QN
div(IQTensor const& T);

IQTensor
combiner(std::vector<IQIndex> inds, Args const& args = Global::args());

IQIndex
combinedIndex(IQTensor const& C);

template<typename... Inds>
IQTensor
combiner(IQIndex const& i1, 
         Inds const&... inds);

//Construct diagonal IQTensor with diagonal 
//elements set to 1.0
template<typename... Inds>
IQTensor
delta(IQIndex const& i1,
      Inds const&... inds);

IQIndex
findIQInd(IQTensor const& T, Index const& i);

QN
qn(IQTensor const& T, Index const& i);

Arrow
dir(IQTensor const& T, Index const& i);

Arrow
dir(IQTensor const& T, IQIndex const& i);

template <typename... IQIndVals>
IQTensor
randomTensor(IQIndexVal const& iv1, 
             IQIndVals&&... ivs);
     
template <typename... Inds>
IQTensor
randomTensor(QN const& q, IQIndex const& i1, Inds &&... inds);

template<typename... Inds>
IQTensor
randomTensor(QN const& q, IQIndexSet const& is);

template <typename... VArgs>
IQTensor
randomTensorC(QN const& q, VArgs&&... vargs);

bool
isEmpty(IQTensor const& T);

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

struct AddITensor;
const char*
typeNameOf(AddITensor const&);

QN inline
flux(IQTensor const& T) { return div(T); }

} //namespace itensor

//See file iqtensor_impl.h for template/inline method implementations
#include "itensor/iqtensor_impl.h"

#endif
