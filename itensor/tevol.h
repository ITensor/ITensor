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
// steps of tstep.
//
// Works by grouping pairs of sites, projecting the MPS
// into the fixed-m manifold, then taking a time step.
//
// Options recognized:
//     Verbose - print useful information to stdout
//
template <class Tensor>
Real
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

template <class Tensor>
void inline
derivMPS(const std::vector<Tensor>& psi, const MPOt<Tensor>& H, 
         std::vector<Tensor>& d, 
         Direction dir = Fromleft)
    {
    DerivMPS<Tensor> D(H,dir);
    d = D(psi);
    }

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
template <class Iterable, class Tensor>
Real
gateTEvol(const Iterable& gatelist, 
          Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts = Global::opts());



//
//
// Implementations
//
//

template <class Iterable, class Tensor>
Real
gateTEvol(const Iterable& gatelist, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts)
    {
    bool verbose = opts.getBool("Verbose",false);

    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    Real tsofar = 0;
    Real tot_norm = psi.normalize();
    if(verbose) 
        {
        Cout << Format("Taking %d steps of timestep %.5f, total time %.5f")
                % nt
                % tstep
                % ttotal
                << Endl;
        }
    for(int tt = 1; tt <= nt; ++tt)
        {
        Foreach(const BondGate<Tensor>& G, gatelist)
            {
            psi.position(G.i());
            psi.applygate(G);
            }

        if(verbose)
            {
            Real percentdone = (100.*tt)/nt;
            if(percentdone < 99.5 || (tt==nt))
                {
                Cout << Format("\b\b\b%2.f%%") % percentdone;
                Cout.flush();
                }
            }

        tot_norm *= psi.normalize();

        tsofar += tstep;
        }
    if(verbose) 
        {
        Cout << Format("\nTotal time evolved = %.5f\n") % tsofar << Endl;
        }

    return tot_norm;

    } // gateTEvol

#undef Cout
#undef Endl
#undef Format

#endif
