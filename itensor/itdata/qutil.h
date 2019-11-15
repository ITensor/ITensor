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

BlOf inline
make_blof(Block const& b, long o)
    {
    BlOf B;
    B.block = b;
    B.offset = o;
    return B;
    }

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

template<typename T>
int
getBlockLoc(QDense<T> const& d,
            Block const& block_ind)
    {
    auto loc = offsetOfLoc(d.offsets,block_ind);
    return loc;
    }

template<typename BlockSparse>
auto
getBlock(BlockSparse & d,
         IndexSet const& is,
         Block const& block_ind)
    -> stdx::if_compiles_return<decltype(makeDataRange(d.data(),d.size())),decltype(d.offsets)>
    {
    auto r = long(block_ind.size());
    if(r == 0) return makeDataRange(d.data(),d.size());
#ifdef DEBUG
    if(is.order() != r) Error("Mismatched size of IndexSet and block_ind in getBlock");
#endif
    //Do binary search to see if there
    //is a block with block index ii
    auto boff = offsetOf(d.offsets,block_ind);
    if(boff >= 0) return makeDataRange(d.data(),boff,d.size());
    using data_range_type = decltype(makeDataRange(d.data(),d.size()));
    return data_range_type{};
    }

// From two input block-sparse tensors,
// output the offsets and data size of the
// result of contracting the tensors
template<typename BlockSparseA,
         typename BlockSparseB>
std::tuple<BlockOffsets,int,std::vector<std::tuple<Block,Block,Block>>>
getContractedOffsets(BlockSparseA const& A,
                     IndexSet const& Ais,
                     BlockSparseB const& B,
                     IndexSet const& Bis,
                     IndexSet const& Cis)
    {
    auto rA = order(Ais);
    auto rB = order(Bis);
    auto rC = order(Cis);

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
    auto Cblockind = Block(rC,0);

    // Store pairs of unordered block numbers and their sizes,
    // to be ordered later
    auto Cblocksizes = BlockOffsets();

    auto blockContractions = std::vector<std::tuple<Block,Block,Block>>();

    //Loop over blocks of A (labeled by elements of A.offsets)
    for(auto const& aio : A.offsets)
        {
        //Begin computing elements of Cblock(=destination of this block-block contraction)
        for(auto iA : range(rA))
            if(AtoC[iA] != -1) Cblockind[AtoC[iA]] = aio.block[iA];

        //Loop over blocks of B which contract with current block of A
        for(auto const& bio : B.offsets)
            {
            auto do_blocks_contract = true;
            for(auto iA : range(rA))
                {
                auto iB = AtoB[iA];
                if(AtoB[iA] != -1)
                    if(aio.block[iA] != bio.block[iB])
                        {
                        do_blocks_contract = false;
                        break;
                        }
                }
            if(!do_blocks_contract) continue;

            //Finish making Cblockind
            for(auto iB : range(rB))
                if(BtoC[iB] != -1) Cblockind[BtoC[iB]] = bio.block[iB];

            // Store the current contraction
            blockContractions.push_back(std::make_tuple(aio.block,bio.block,Cblockind));

            long blockDim = 1;   //accumulate dim of Indices
            for(auto j : range(order(Cis)))
                {
                auto& J = Cis[j];
                auto i_j = Cblockind[j];
                blockDim *= J.blocksize0(i_j);
                }

            Cblocksizes.push_back(make_blof(Cblockind,blockDim));
            } //for B.offsets
        } //for A.offsets

    // Sort the block sizes by the block labels
    std::sort(Cblocksizes.begin(),Cblocksizes.end(),
              [](auto a, auto b) { return a.block < b.block; });

    // Remove the duplicates, need to resize manually
    auto newCend = std::unique(Cblocksizes.begin(),
                               Cblocksizes.end(),
                               [](auto a, auto b) { return a.block == b.block; } );
    Cblocksizes.resize(std::distance(Cblocksizes.begin(),newCend));

    auto current_offset = 0;
    for(auto i : range(Cblocksizes.size()))
        {
        auto current_size = Cblocksizes[i].offset;
        Cblocksizes[i].offset = current_offset;
        current_offset += current_size;
        } 
    // Stores the total size that the storage of C should have
    auto Csize = current_offset;

    return std::make_tuple(Cblocksizes,Csize,blockContractions);
    }

template<typename TA,
         typename TB,
         typename TC,
         typename Callable>
void
loopContractedBlocks(QDense<TA> const& A,
                     IndexSet const& Ais,
                     QDense<TB> const& B,
                     IndexSet const& Bis,
                     QDense<TC> & C,
                     IndexSet const& Cis,
                     std::vector<std::tuple<Block,Block,Block>> const& blockContractions,
                     Callable & callback)
    {
    for(auto const& [Ablockind,Bblockind,Cblockind] : blockContractions)
        {
        auto ablock = getBlock(A,Ais,Ablockind);
        auto bblock = getBlock(B,Bis,Bblockind);
        auto cblock = getBlock(C,Cis,Cblockind);
        auto Cblockloc = getBlockLoc(C,Cblockind);
        callback(ablock,Ablockind,
                 bblock,Bblockind,
                 cblock,Cblockind,
                 Cblockloc);
        }
    }

// This is a special case of loopContractedBlocks for QDiag
// since QDiag doesn't have a .offsets function
// TODO: add .offsets for QDiag so it can use the generic
// (faster) version
template<typename TA, 
         typename TB,
         typename TC,
         typename Callable>
void
loopContractedBlocks(TA const& A,
                     IndexSet const& Ais,
                     TB const& B,
                     IndexSet const& Bis,
                     TC & C,
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
    auto Bblockind = Block(rB,0);
    auto Cblockind = Block(rC,0);
    //Loop over blocks of A (labeled by elements of A.offsets)
    for(auto const& aio : A.offsets)
        {
        //Reset couB to run over indices of B (at first)
        couB.reset();
        for(auto iB : range(rB))
            {
            couB.setRange(iB,0,Bis[iB].nblock()-1);
            }
        for(auto iA : range(rA))
            {
            auto ival = aio.block[iA];
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
            //Finish making Cblockind and Bblockind
            for(auto iB : range(rB))
                {
                if(BtoC[iB] != -1) Cblockind[BtoC[iB]] = couB.i[iB];
                Bblockind[iB] = couB.i[iB];
                }
            auto bblock = getBlock(B,Bis,Bblockind);
            if(!bblock) continue;

            auto cblock = getBlock(C,Cis,Cblockind);
            assert(cblock);

            auto ablock = makeDataRange(A.data(),aio.offset,A.size());

            callback(ablock,aio.block,
                     bblock,Bblockind,
                     cblock,Cblockind);
            } //for couB
        } //for A.offsets
    }


} //namespace itensor

#endif
