//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TDVP_H
#define __ITENSOR_TDVP_H

#include "mpo.h"

template <class Tensor>
void
imagTEvol(const MPOt<Tensor>& H, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi);

template <class Tensor>
void
exactImagTEvol(const MPOt<Tensor>& H, Real ttotal, Real tstep, 
               MPSt<Tensor>& psi);

#endif
