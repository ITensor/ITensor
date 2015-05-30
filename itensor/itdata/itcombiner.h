//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITCOMBINER_H
#define __ITENSOR_ITCOMBINER_H

#include "itensor/itdata/itreal.h"

namespace itensor {

class ITCombiner : public RegisterData<ITCombiner>
    {
    public:

    ITCombiner() { }

    virtual
    ~ITCombiner() { }

    };

void inline
read(std::istream& s, ITCombiner& dat) { }

void inline
write(std::ostream& s, const ITCombiner& dat) { }

Cplx
doTask(const GetElt<Index>& g, const ITCombiner& c);

Real
doTask(const NormNoScale<Index>& N, const ITCombiner& d);

void
doTask(Conj,const ITCombiner& d);

void
doTask(Contract<Index>& C,
       const ITReal& d,
       const ITCombiner& cmb,
       ManagePtr& mp);

void
doTask(Contract<Index>& C,
       const ITCombiner& cmb,
       const ITReal& d,
       ManagePtr& mp);

bool
doTask(CheckComplex, const ITCombiner& d);

void
doTask(PrintIT<Index>& P, const ITCombiner& d);

void
doTask(Write& W, const ITCombiner& d);

} //namespace itensor

#endif

