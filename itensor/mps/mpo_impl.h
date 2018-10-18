//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPO_IMPL_H
#define __ITENSOR_MPO_IMPL_H

namespace itensor {

template <> inline
MPO IQMPO::
toMPO() const
    {
    auto res = MPO(sites_,logrefNorm_);
    for(int j = 0; j <= N()+1; ++j)
        {
        //if(!res.A_.at(j)) continue;
        res.A_.at(j) = toITensor(A(j));
        }
    return res;
    }
 
//toMPO method fails unless template class 
//Tensor is set to IQTensor (object is an IQMPO)
template<class Tensor>
MPO MPOt<Tensor>::
toMPO() const
    {
    Error("toMPO only implemented for class IQMPO");
    return MPO();
    }


template<typename T>
bool
isOrtho(MPOt<T> const& W)
    {
    return W.leftLim()+1 == W.rightLim()-1;
    }

template<typename T>
int
orthoCenter(MPOt<T> const& W)
    {
    if(!isOrtho(W)) Error("orthogonality center not well defined.");
    return (W.leftLim() + 1);
    }

template <class Tensor>
MPOt<Tensor>
sum(MPOt<Tensor> L, 
    MPOt<Tensor> const& R, 
    Args const& args)
    {
    L.plusEq(R,args);
    return L;
    }

template <class Tensor>
void 
psiHphi(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        MPSt<Tensor> const& phi, 
        Real& re, 
        Real& im)
    {
    overlap(psi,H,phi,re,im);
    }

template <class Tensor>
Real 
psiHphi(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        MPSt<Tensor> const& phi) //Re[<psi|H|phi>]
    {
    return overlap(psi,H,phi);
    }

template <class Tensor>
Complex 
psiHphiC(MPSt<Tensor> const& psi, 
         MPOt<Tensor> const& H, 
         MPSt<Tensor> const& phi) //Re[<psi|H|phi>]
    {
    return overlapC(psi,H,phi);
    }

template<class Tensor>
void
psiHphi(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        Tensor const& LB, 
        Tensor const& RB, 
        MPSt<Tensor> const& phi, 
        Real& re, 
        Real& im) //<psi|H|phi>
    {
    overlap(psi,H,LB,RB,phi,re,im);
    }

template <class Tensor>
Real
psiHphi(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        Tensor const& LB, 
        Tensor const& RB, 
        MPSt<Tensor> const& phi) //Re[<psi|H|phi>]
    {
    return overlap(psi,H,LB,RB,phi);
    }

template <class Tensor>
void
psiHKphi(MPSt<Tensor> const& psi, 
         MPOt<Tensor> const& H, 
         MPOt<Tensor> const& K,
         MPSt<Tensor> const& phi, 
         Real& re, 
         Real& im) //<psi|H K|phi>
    {
    overlap(psi,H,K,phi,re,im);
    }

template <class Tensor>
Real
psiHKphi(MPSt<Tensor> const& psi, 
         MPOt<Tensor> const& H, 
         MPOt<Tensor> const& K,
         MPSt<Tensor> const& phi) //<psi|H K|phi>
    {
    return overlap(psi,H,K,phi);
    }

template <class Tensor>
Complex
psiHKphiC(MPSt<Tensor> const& psi, 
          MPOt<Tensor> const& H, 
          MPOt<Tensor> const& K,
          MPSt<Tensor> const& phi) //<psi|H K|phi>
    {
    return overlapC(psi,H,K,phi);
    }

} //namespace itensor

#endif
