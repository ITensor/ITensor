//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEVOL_H
#define __ITENSOR_TEVOL_H

#include "mpo.h"
#include "bondgate.h"
#include <list>

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

//
// Imaginary time evolve an MPS by an amount ttotal in time
// steps of tstep using the Hamiltonian MPO H.
//
// Options recognized:
//     Verbose - print useful information to stdout
//     Order - order at which to stop applying powers of H, 
//             setting order to p yields error of tstep^p
//     Maxm - Maximum states kept each step
//     Cutoff - Maximum truncation error each step
//     Nsweep - Number of sweeps used to apply H to MPS (see fitApplyMPO)
//
template <class Tensor>
void
imagTEvol(const MPOt<Tensor>& H, 
          Real ttotal, 
          Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts = Global::opts());


//
// Evolves an MPS in real or imaginary time by an amount ttotal in steps
// of tstep using the list of bond gates provided.
//
// Options recognized:
//     Verbose - print useful information to stdout
//
template <class Tensor>
Real
gateTEvol(const std::list<BondGate<Tensor> >& gatelist, 
          Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts = Global::opts());


#undef Cout
#undef Endl
#undef Format

#endif
