//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BONDGATE_H
#define __ITENSOR_BONDGATE_H

#include "iqtensor.h"

template <class Tensor>
class BondGate
    {
    public:

    enum Type { tReal,  //real-time gate
                tImag,  //imaginary-time gate
                Swap }; //exchange states of site i and j

    BondGate(const Model& model, int i, int j);

    BondGate(const Model& model, int i, int j, 
             Type type, Real tau, Tensor bondH);

    operator const Tensor&() const { return gate_; }

    const Tensor&
    gate() const { return gate_; }

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
BondGate<Tensor>::
BondGate(const Model& model, int i, int j)
    : 
    type_(Swap), 
    i_(i), 
    j_(j)
    {
    Tensor a(model.si(i),primed(model.si(j)),1);
    Tensor b(model.si(j),primed(model.si(i)),1);
    gate_ = a*b;
    }

template <class Tensor>
BondGate<Tensor>::
BondGate(const Model& model, int i, int j, 
         Type type, Real tau, Tensor bondH)
    : 
    type_(type), 
    i_(i), 
    j_(j)
    {
    if(!(type_ == tReal || type_ ==tImag))
        {
        Error("When providing bondH, type must be tReal or tImag");
        }
    bondH *= -tau;
    Tensor unit = model.id(i)*model.id(j);
    if(type_ == tReal)
        {
        bondH *= Complex_i;
        }
    Tensor term = bondH;
    bondH.mapprime(1,2);
    bondH.mapprime(0,1);

    // exp(x) = 1 + x +  x^2/2! + x^3/3! ..
    // = 1 + x * (1 + x/2 *(1 + x/3 * (...
    // ~ ((x/3 + 1) * x/2 + 1) * x + 1
    for(int ord = 100; ord >= 1; --ord)
        {
        term /= ord;
        gate_ = unit + term;
        term = gate_ * bondH;
        term.mapprime(2,1);
        }
    }

#endif
