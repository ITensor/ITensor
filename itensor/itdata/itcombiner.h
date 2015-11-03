//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITCOMBINER_H
#define __ITENSOR_ITCOMBINER_H

#include "itensor/itdata/dense.h"

namespace itensor {

class QN;

class ITCombiner
    {
    public:

    ITCombiner() { }

    };

void inline
read(std::istream& s, ITCombiner& dat) { }

void inline
write(std::ostream& s, const ITCombiner& dat) { }

Cplx
doTask(const GetElt<Index>& g, const ITCombiner& c);

Real
doTask(NormNoScale, const ITCombiner& d);

void
doTask(Conj,const ITCombiner& d);

template<typename V>
void
doTask(Contract<Index> & C,
       Dense<V>   const& d,
       ITCombiner const& cmb,
       ManageStore     & m);

template<typename V>
void
doTask(Contract<Index> & C,
       ITCombiner const& cmb,
       Dense<V>   const& d,
       ManageStore     & m);

bool
doTask(CheckComplex, ITCombiner const& d);

void
doTask(PrintIT<Index>& P, ITCombiner const& d);

auto inline
doTask(StorageType const& S, ITCombiner const& d) ->StorageType::Type { return StorageType::ITCombiner; }

QN 
doTask(CalcDiv const& C, ITCombiner const& d);

} //namespace itensor

#endif

