//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_FERMION_H
#define __ITENSOR_FERMION_H
#include "itensor/iqindex.h"
#include "itensor/tensor/permutation.h"

namespace itensor {

// from two indexsets (possibly) sharing some indices, produce the permutations
// of their labels (0-indexed) which would produce 'standard' orientation e.g.
// if:      Ais   = J,K,I               Bis   = P,I,J,Q
// then:    permA = (0,1),(1,0),(2,2)   permB = (0,2),(1,1),(2,0),(3,3)
template <class index_type>
void inline
unriffle(IndexSetT<index_type> const& Ais,
         IndexSetT<index_type> const& Bis,
         Permutation & permA,
         Permutation & permB)
    {
    auto rA = Ais.r();
    auto rB = Bis.r();
    int   n = 0;
    int   k = 0;
    permA   = Permutation(rA,-1);   //reinitialize the permutations
    permB   = Permutation(rB,-1);   //with default values -1

    for(auto i:range(rA))           //ascending i= 0,...,rA-1
        {
        for(auto j: range(rB))      //ascending j=0,...,rB-1
            {
            if(Ais[i]==Bis[j])
                {
                permA[i] = n;
                permB[j] = n;
                n++;
                continue;
                }
            }
        }

    for(auto i: range(rA))
        {
        if(permA.dest(i)==-1)
            {
            permA[i]=k;
            k++;
            }
        else
            {
            permA[i]=permA[i]+rA-n;
            }
        }

    for(auto j:range(rB))
        {
        if(permB.dest(j)==-1)
            {
            permB[j]=n;
            n++;
            }
        }
    }
template void unriffle(IQIndexSet const& Ais,
                       IQIndexSet const& Bis,
                       Permutation & permA,
                       Permutation & permB);


// Given a permutation P and an ordered subset of P's domain R,
// returns the number of inversions in a sort of the image
// of P precomposed with R:
//      r0 r1 ... rn ---> P[r0] P[r1] ... P[rn]
// When R is ordered, returns number of crossings in P restricted
// to the domain R.
// mod 2 this is the signature of P|R.
int inline
count_filtered_swaps(Permutation const& P, IntArray R)	//ask Miles, why fail for &R ?
    {
    auto rR = R.size();
    for(auto i:range(rR))
        {
        R[i]=P[R[i]];
        }
	
	int swaps = 0;
    auto Y = IntArray(P.size(),0);
    for(int i=rR-1; i>=0; --i)
        {
        auto x = R[i];
        for(auto j: range(x))
            {
            if(Y[j] != 0){swaps++;}
            }
        Y[x]=1;
        }

    return swaps;
    }

int inline
total_swaps(Permutation const& P, IQIndexSet const I, IntArray const& block)
    {
    //probably add an assert I.size() == block.size() == P.size()
    auto rI = I.r();
    int num_swaps = 0;

    //assumes ALL IQIndex have same total space of QN's
    //use the first index to determine to determine where different
    //'colors' of fermions can be.
    auto i1 = I[0].qn(1);  //first index/sector of 0th IQIndex of B
    for(auto c=1; c<=QNSize(); ++c)
        {
        if(not isActive(i1,c)) break;        //returns 1 until past dim of qn
        if(not isFermionic(i1,c)) continue;  //skip non-fermionic
        auto cib = IntArray();
        for(auto i:range(rI))
            {
            //value of the quantum number of the c colored fermion
            //within this block of the tensor
            auto q = I[i].qn(block[i]+1)[c-1];
            if(std::abs(q)%2==1)
                {
                cib.push_back(i);
                }
            }
        num_swaps += count_filtered_swaps(P,cib);
        }
    return num_swaps;
    }


} //namespace itensor

#endif