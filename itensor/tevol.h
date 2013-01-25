//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TDVP_H
#define __ITENSOR_TDVP_H

#include "mpo.h"

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


//Forward declaration
template <class Tensor> class BondGate;

//
// Evolves an MPS in real or imaginary time by an amount ttotal in steps
// of tstep using the list of bond gates provided.
//
// Options recognized:
//     Verbose - print useful information to stdout
//
template <class Tensor>
void
gateTEvol(const std::list<BondGate<Tensor> >& gatelist, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts = Global::opts());

template <class Tensor>
class BondGate
    {
    public:

    enum Type { ExpH, SwapOp };

    BondGate(const Model& model, int i, int j)
        : 
        type_(SwapOp), 
        i_(i), 
        j_(j)
        {
        Tensor a(model.si(i),primed(model.si(j)),1);
        Tensor b(model.si(j),primed(model.si(i)),1);
        gate_ = a*b;
        }

    BondGate(const Model& model,Type t, int i, int j, 
             Real tau, Tensor bondH)
        : 
        type_(ExpH), 
        i_(i), 
        j_(j)
        {
        bondH *= -tau;
        Tensor unit = model.id(i)*model.id(j);
        Tensor term = bondH, kk;
        bondH.mapprime(1,2);
        bondH.mapprime(0,1);

        // exp(x) = 1 + x +  x^2/2! + x^3/3! ..
        // = 1 + x * (1 + x/2 *(1 + x/3 * (...
        // ~ ((x/3 + 1) * x/2 + 1) * x + 1
        for(int o = 100; o >= 1; o--)
            {
            if(o != 1) term *= 1.0 / o;
            kk = unit + term;
            term = kk * bondH;
            term.mapprime(2,1);
            }
        gate_ = kk;
        }

    const Tensor&
    op() const { return gate_; }

    int
    i() const { return i_; }

    int
    j() const { return j_; }

    Type
    type() const { return type_; }

    private:

    Type type_;
    int i_,j_; // The left, right indices of bond
    Tensor gate_;

    };

template <class Tensor>
void
derivMPS(const std::vector<Tensor>& psi, const MPOt<Tensor>& H, 
         std::vector<Tensor>& dpsi, Direction dir = Fromleft);

#endif
