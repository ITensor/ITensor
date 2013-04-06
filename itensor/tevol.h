//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEVOL_H
#define __ITENSOR_TEVOL_H

#include "mpo.h"
#include "bondgate.h"

//
// Imaginary time evolve an MPS by an amount ttotal in time
// steps of tstep.
//
// Works by grouping pairs of sites, projecting the MPS
// into the fixed-m manifold, then taking a time step.
//
// Options recognized:
//     Verbose - print useful information to stdout
//
template <class Tensor>
void
imagTEvol(const MPOt<Tensor>& H, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts = Global::opts());

template <class Tensor>
void
derivMPS(const std::vector<Tensor>& psi, const MPOt<Tensor>& H, 
         std::vector<Tensor>& dpsi, Direction dir = Fromleft);


//
// Evolves an MPS in real or imaginary time by an amount ttotal in steps
// of tstep using the list of bond gates provided.
//
// Options recognized:
//     Verbose - print useful information to stdout
//
template <class Tensor>
void
gateTEvol(const std::list<BondGate<Tensor> >& gatelist, 
          Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts = Global::opts());

#endif
