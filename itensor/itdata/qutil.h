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
#ifndef __ITENSOR_QUTIL_H
#define __ITENSOR_QUTIL_H

#include "itensor/indexset.h"

namespace itensor {


//
// Helper object for treating
// QDense storage as a "tensor of tensors"
//
template<typename Indexable>
class IndexDim
    {
    IndexSet const& is_;
    Indexable const& ind_;
    public:

    IndexDim(IndexSet const& is,
             Indexable const& ind)
      : is_(is),
        ind_(ind)
        { }

    size_t
    size() const { return is_.order(); }

    //size_t
    //operator[](size_t j) const { return (is_[j])[dim(ind_[j]]); }
    size_t
    operator[](size_t j) const { return (is_[j]).blocksize0(ind_[j]); }
    };

template<typename Indexable>
auto
make_indexdim(IndexSet const& is, Indexable const& ind) 
    -> IndexDim<Indexable>
    { 
    return IndexDim<Indexable>(is,ind); 
    }

// For a block index (0,1,...,Nblocks-1),
// as in the offsets array of an QDense,
// computes the zero-indexed location of
// the given block (e.g. 5 -> {1,0,2})
// storing the resulting indices in the 
// container "ind". The jth index of ind
// can range from 0 to is[j].nblock(), 
// such that these indices correspond to
// viewing the QDense storage as a 
// "tensor of tensors"
template<typename Container>
void
computeBlockInd(long block,
                IndexSet const& is,
                Container & ind)
    {
    using size_type = decltype(ind.size());
    size_type r = ind.size();
    assert(r == size_type(is.order()));
    for(size_type j = 0; j < r-1; ++j)
        {
        auto res = std::ldiv(block,is[j].nblock());
        ind[j] = res.rem;
        block = res.quot;
        }
    ind[r-1] = block;
    }

template<typename BlockSparse, typename Indexable>
auto
getBlock(BlockSparse & d,
         IndexSet const& is,
         Indexable const& block_ind)
    -> stdx::if_compiles_return<decltype(makeDataRange(d.data(),d.size())),decltype(d.offsets)>
    {
    auto r = long(block_ind.size());
    if(r == 0) return makeDataRange(d.data(),d.size());
#ifdef DEBUG
    if(is.order() != r) Error("Mismatched size of IndexSet and block_ind in getBlock");
#endif
    long ii = 0;
    for(auto i = r-1; i > 0; --i)
        {
        ii += block_ind[i];
        ii *= is[i-1].nblock();
        }
    ii += block_ind[0];
    //Do binary search to see if there
    //is a block with block index ii
    auto boff = offsetOf(d.offsets,ii);
    if(boff >= 0) return makeDataRange(d.data(),boff,d.size());
    using data_range_type = decltype(makeDataRange(d.data(),d.size()));
    return data_range_type{};
    }

template<typename BlockSparseA, 
         typename BlockSparseB,
         typename BlockSparseC,
         typename Callable>
void
loopContractedBlocks(BlockSparseA const& A,
                     IndexSet const& Ais,
                     BlockSparseB const& B,
                     IndexSet const& Bis,
                     BlockSparseC & C,
                     IndexSet const& Cis,
                     Callable & callback)
    {
    auto rA = Ais.order();
    auto rB = Bis.order();
    auto rC = Cis.order();

    auto AtoB = IntArray(rA,-1);
    auto AtoC = IntArray(rA,-1);
    auto BtoC = IntArray(rB,-1);
    for(auto ic : range(rC))
        {
        auto j = indexPosition(Ais,Cis[ic]);
        if(j >= 0)
            {
            AtoC[j] = ic;
            }
        else
            {
            j = indexPosition(Bis,Cis[ic]);
            BtoC[j] = ic;
            }
        }
    for(auto ia : range(rA))
    for(auto ib : range(rB))
        {
        if(Ais[ia] == Bis[ib])
            {
            AtoB[ia] = ib;
            break;
            }
        }

    auto couB = detail::GCounter(rB);
    auto Ablockind = IntArray(rA,0);
    auto Cblockind = IntArray(rC,0);
    //Loop over blocks of A (labeled by elements of A.offsets)
    for(auto& aio : A.offsets)
        {
        //Reconstruct indices labeling this block of A, put into Ablock
        //TODO: optimize away need to call computeBlockInd by
        //      storing block indices directly in QDense
        //      Taking 10% of running time in S=1 N=100 DMRG tests (maxdim=100)
        computeBlockInd(aio.block,Ais,Ablockind);
        //Reset couB to run over indices of B (at first)
        couB.reset();
        for(auto iB : range(rB))
            {
            couB.setRange(iB,0,Bis[iB].nblock()-1);
            }
        for(auto iA : range(rA))
            {
            auto ival = Ablockind[iA];
            //Restrict couB to be fixed for indices of B contracted with A
            if(AtoB[iA] != -1) couB.setRange(AtoB[iA],ival,ival);
            //Begin computing elements of Cblock(=destination of this block-block contraction)
            if(AtoC[iA] != -1) Cblockind[AtoC[iA]] = ival;
            }
        //Loop over blocks of B which contract with current block of A
        for(;couB.notDone(); ++couB)
            {
            //Check whether B contains non-zero block for this setting of couB
            //TODO: check whether block is present by storing all blocks
            //      but most have null pointers to data
            auto bblock = getBlock(B,Bis,couB.i);
            if(!bblock) continue;

            //Finish making Cblockind and Bblockind
            auto Bblockind = IntArray(rB,0);
            for(auto iB : range(rB))
                {
                if(BtoC[iB] != -1) Cblockind[BtoC[iB]] = couB.i[iB];
                Bblockind[iB] = couB.i[iB];
                }

            auto cblock = getBlock(C,Cis,Cblockind);
            assert(cblock);

            auto ablock = makeDataRange(A.data(),aio.offset,A.size());

            callback(ablock,Ablockind,
                     bblock,Bblockind,
                     cblock,Cblockind);

            } //for couB
        } //for A.offsets
    }


} //namespace itensor

#endif
