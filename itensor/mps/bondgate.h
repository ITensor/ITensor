//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BONDGATE_H
#define __ITENSOR_BONDGATE_H

#include "itensor/iqtensor.h"
#include "itensor/mps/siteset.h"

namespace itensor {

template <class Tensor>
class BondGate;

using Gate = BondGate<ITensor>;
using IQGate = BondGate<IQTensor>;

template <class Tensor>
class BondGate
    {
    public:

    enum Type { tReal,  //real-time gate
                tImag,  //imaginary-time gate
                Swap,
                Custom }; //exchange states of sites i1 and i2

    BondGate(SiteSet const& sites, 
             int i1, 
             int i2);

    BondGate(SiteSet const& sites, 
             int i1, 
             int i2, 
             Type type, 
             Real tau, 
             Tensor bondH);

    BondGate(SiteSet const& sites, 
             int i1, 
             int i2,
             Tensor gate);

    int i1() const { return i1_; }

    int i2() const { return i2_; }

    operator const Tensor&() const { return gate_; }

    Tensor const&
    gate() const { return gate_; }

    Type
    type() const { return type_; }

    private:

    Type type_;
    int i1_,i2_; // The left, right indices of bond
    Tensor gate_;

    void
    makeSwapGate(SiteSet const& sites);
    };

template<class Tensor>
Tensor
operator*(BondGate<Tensor> const& G, Tensor T) { T *= G.gate(); return T; }

template<class Tensor>
Tensor
operator*(Tensor T, BondGate<Tensor> const& G) { T *= G.gate(); return T; }

template <class Tensor>
BondGate<Tensor>::
BondGate(SiteSet const& sites, 
         int i1, 
         int i2)
  : type_(Swap) 
    {
    if(i1 < i2)
        {
        i1_ = i1;
        i2_ = i2;
        }
    else
        {
        i1_ = i2;
        i2_ = i1;
        }
    makeSwapGate(sites);
    }

template <class Tensor>
BondGate<Tensor>::
BondGate(SiteSet const& sites, 
         int i1, 
         int i2, 
         Type type, 
         Real tau, 
         Tensor bondH)
  : type_(type)
    {
    if(i1 < i2)
        {
        i1_ = i1;
        i2_ = i2;
        }
    else
        {
        i1_ = i2;
        i2_ = i1;
        }

    if(!(type_ == tReal || type_ ==tImag))
        {
        Error("When providing bondH, type must be tReal or tImag");
        }
    bondH *= -tau;
    Tensor unit = sites.op("Id",i1_)*sites.op("Id",i2_);
    if(type_ == tReal)
        {
        bondH *= Complex_i;
        }
    auto term = bondH;
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
BondGate<Tensor>::
BondGate(SiteSet const& sites, 
         int i1, 
         int i2,
         Tensor gate)
  : type_(Custom)
    {
    i1_ = i1;
    i2_ = i2;
    gate_ = gate;
    }

template <class Tensor>
void BondGate<Tensor>::
makeSwapGate(SiteSet const& sites)
    {
    auto s1 = sites(i1_);
    auto s2 = sites(i2_);
    auto a = Tensor(dag(s1),prime(s2));
    auto b = Tensor(dag(s2),prime(s1));
    for(auto j : range1(s1))
        {
        a.set(dag(s1)(j),prime(s2)(j),1.);
        b.set(dag(s2)(j),prime(s1)(j),1.);
        }
    gate_ = a*b;
    }

} //namespace itensor

#endif
