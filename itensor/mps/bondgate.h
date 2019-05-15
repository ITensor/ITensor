//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_BONDGATE_H
#define __ITENSOR_BONDGATE_H

#include "itensor/itensor.h"
#include "itensor/mps/siteset.h"

namespace itensor {

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
             ITensor bondH);

    BondGate(SiteSet const& sites, 
             int i1, 
             int i2,
             ITensor gate);

    int i1() const { return i1_; }

    int i2() const { return i2_; }

    operator const ITensor&() const { return gate_; }

    ITensor const&
    gate() const { return gate_; }

    Type
    type() const { return type_; }

    private:

    Type type_;
    int i1_,i2_; // The left, right indices of bond
    ITensor gate_;

    void
    makeSwapGate(SiteSet const& sites);
    };

ITensor inline
operator*(BondGate const& G, ITensor T) { T *= G.gate(); return T; }

ITensor inline
operator*(ITensor T, BondGate const& G) { T *= G.gate(); return T; }

inline BondGate::
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

inline BondGate::
BondGate(SiteSet const& sites, 
         int i1, 
         int i2, 
         Type type, 
         Real tau, 
         ITensor bondH)
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
    ITensor unit = sites.op("Id",i1_)*sites.op("Id",i2_);
    if(type_ == tReal)
        {
        bondH *= Complex_i;
        }
    auto term = bondH;
    bondH.replaceTags("1","2");
    bondH.replaceTags("0","1");

    // exp(x) = 1 + x +  x^2/2! + x^3/3! ..
    // = 1 + x * (1 + x/2 *(1 + x/3 * (...
    // ~ ((x/3 + 1) * x/2 + 1) * x + 1
    for(int ord = 100; ord >= 1; --ord)
        {
        term /= ord;
        gate_ = unit + term;
        term = gate_ * bondH;
        term.replaceTags("2","1");
        }
    }

inline BondGate::
BondGate(SiteSet const& sites, 
         int i1, 
         int i2,
         ITensor gate)
  : type_(Custom)
    {
    i1_ = i1;
    i2_ = i2;
    gate_ = gate;
    }

void inline BondGate::
makeSwapGate(SiteSet const& sites)
    {
    auto s1 = sites(i1_);
    auto s2 = sites(i2_);
    auto a = ITensor(dag(s1),prime(s2));
    auto b = ITensor(dag(s2),prime(s1));
    for(auto j : range1(s1))
        {
        a.set(dag(s1)(j),prime(s2)(j),1.);
        b.set(dag(s2)(j),prime(s1)(j),1.);
        }
    gate_ = a*b;
    }

} //namespace itensor

#endif
