//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEVOL_H
#define __ITENSOR_TEVOL_H

#include "mpo.h"
#include "bondgate.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

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
class DerivMPS
    {
    public:

    DerivMPS(const MPOt<Tensor>& H, Direction dir = Fromleft)
        : H_(H), dir_(dir) { }
    
    std::vector<Tensor>
    operator()(const std::vector<Tensor>& psi) const;

    private:

    const MPOt<Tensor>& H_;
    const Direction dir_;

    };

//
// Compute the norm (= sqrt(|<psi|psi>|)) of an
// MPS-like vector of tensors
// 
// vector psi is 1-indexed
// Automatically determines size by counting number
// of non-Null tensors
//
template <class Tensor>
Real
norm(const std::vector<Tensor>& psi);

//
// Compute the expectation value (= <psi|H|psi>)
// of an MPS-like vector of tensors with respect 
// to an MPO H
//
// vector psi is 1-indexed
// Automatically determines size by counting number
// of non-Null tensors
//
template <class Tensor>
Real
expect(const std::vector<Tensor>& psi, const MPOt<Tensor>& H);

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


#undef Cout
#undef Endl
#undef Format

#endif
