//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BONDGATE_H
#define __ITENSOR_BONDGATE_H

#include "iqtensor.h"
#include "model.h"

template <class Tensor>
class BondGate;

typedef BondGate<ITensor>
Gate;

typedef BondGate<IQTensor>
IQGate;

template <class Tensor>
class BondGate
    {
    public:

    enum Type { tReal,  //real-time gate
                tImag,  //imaginary-time gate
                Swap }; //exchange states of site i and j

    BondGate(const Model& sites, int i, int j);

    BondGate(const Model& sites, int i, int j, 
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

    void
    makeSwapGate(const Model& sites);

    };

template <class Tensor>
BondGate<Tensor>::
BondGate(const Model& sites, int i, int j)
    : 
    type_(Swap), 
    i_(i), 
    j_(j)
    {
    makeSwapGate(sites);
    }

template <class Tensor>
BondGate<Tensor>::
BondGate(const Model& sites, int i, int j, 
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
    Tensor unit = sites.op("Id",i)*sites.op("Id",j);
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

template <class Tensor>
void BondGate<Tensor>::
makeSwapGate(const Model& sites)
    {
    Tensor a(sites(i_),primed(sites(j_)),1);
    Tensor b(sites(j_),primed(sites(i_)),1);
    gate_ = a*b;
    }

template<>
void BondGate<IQTensor>::
makeSwapGate(const Model& sites)
    {
    IQTensor a(conj(sites(i_)),primed(sites(j_))),
             b(conj(sites(j_)),primed(sites(i_)));
    for(int n = 1; n <= sites(i_).nindex(); ++n)
        {
        const Index &iind(sites(i_).index(n)),
                    &jind(sites(j_).index(n));
        a += ITensor(iind,primed(jind),1);
        b += ITensor(jind,primed(iind),1);
        }
    gate_ = a*b;
    }

#endif
