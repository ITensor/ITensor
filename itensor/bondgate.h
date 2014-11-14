//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BONDGATE_H
#define __ITENSOR_BONDGATE_H

#include "iqtensor.h"
#include "siteset.h"

namespace itensor {

template <class Tensor>
class BondGate
    {
    public:

    enum Type { tReal,  //real-time gate
                tImag,  //imaginary-time gate
                Swap }; //exchange states of sites i1 and i2

    BondGate(const Model& sites, int i1, int i2);

    BondGate(const Model& sites, int i1, int i2, 
             Type type, Real tau, Tensor bondH);

    int i1() const { return i1_; }

    int i2() const { return i2_; }

    operator const Tensor&() const { return gate_; }

    const Tensor&
    gate() const { return gate_; }

    Type
    type() const { return type_; }

  



    // Deprecated: use i1() and i2() instead
    int
    i() const { return i1_; }

    // Deprecated: use i1() and i2() instead
    int
    j() const { return i2_; }


    private:

    Type type_;
    int i1_,i2_; // The left, right indices of bond
    Tensor gate_;

    void
    makeSwapGate(const Model& sites);


    };
using Gate = BondGate<ITensor>;
using IQGate = BondGate<IQTensor>;

template<class Tensor>
Tensor
operator*(const BondGate<Tensor>& G, Tensor T) { T *= G.gate(); return T; }

template<class Tensor>
Tensor
operator*(Tensor T, const BondGate<Tensor>& G) { T *= G.gate(); return T; }

template <class Tensor>
BondGate<Tensor>::
BondGate(const Model& sites, int i1, int i2)
    : 
    type_(Swap) 
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
BondGate(const Model& sites, int i1, int i2, 
         Type type, Real tau, Tensor bondH)
    : 
    type_(type)
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
    Tensor a(sites(i1_),prime(sites(i2_)),1);
    Tensor b(sites(i2_),prime(sites(i1_)),1);
    gate_ = a*b;
    }

template<>
void inline BondGate<IQTensor>::
makeSwapGate(const Model& sites)
    {
    IQTensor a(dag(sites(i1_)),prime(sites(i2_))),
             b(dag(sites(i2_)),prime(sites(i1_)));
    for(int n = 1; n <= sites(i1_).nindex(); ++n)
        {
        const Index &i1ind(sites(i1_).index(n)),
                    &i2ind(sites(i2_).index(n));
        a += ITensor(i1ind,prime(i2ind),1);
        b += ITensor(i2ind,prime(i1ind),1);
        }
    gate_ = a*b;
    }

}; //namespace itensor

#endif
