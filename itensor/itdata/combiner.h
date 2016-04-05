//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITCOMBINER_H
#define __ITENSOR_ITCOMBINER_H

#include "itensor/itdata/dense.h"

namespace itensor {

class QN;

class Combiner
    {
    public:

    Combiner() { }

    };

const char*
typeNameOf(Combiner const& d);

void inline
read(std::istream& s, Combiner& dat) { }

void inline
write(std::ostream& s, Combiner const& dat) { }

Cplx
doTask(GetElt<Index> const& g, Combiner const& c);

Real
doTask(NormNoScale, Combiner const& d);

void
doTask(Conj, Combiner const& d);

template<typename V>
void
doTask(Contract<Index> & C,
       Dense<V>   const& d,
       Combiner const& cmb,
       ManageStore     & m);

template<typename V>
void
doTask(Contract<Index> & C,
       Combiner const& cmb,
       Dense<V>   const& d,
       ManageStore     & m);

bool
doTask(CheckComplex, Combiner const& d);

void
doTask(PrintIT<Index>& P, Combiner const& d);

auto inline
doTask(StorageType const& S, Combiner const& d) 
    ->StorageType::Type 
    { return StorageType::Combiner; }

QN 
doTask(CalcDiv const& C, Combiner const& d);

} //namespace itensor

#endif

