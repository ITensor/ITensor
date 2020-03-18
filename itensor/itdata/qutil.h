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

#ifdef ITENSOR_USE_OMP
#include <omp.h>
#endif

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

    // Store pairs of unordered block numbers and their sizes,
    // to be ordered later
    using BlockContractions = std::vector<std::tuple<Block,Block,Block>>;

    //Loop over blocks of A (labeled by elements of A.offsets)
#ifdef ITENSOR_USE_OMP
    int num_threads = omp_get_max_threads();
    auto Cblocksizes_thread = std::vector<BlockOffsets>(num_threads);
    auto blockContractions_thread = std::vector<BlockContractions>(num_threads);
#else
    auto Cblocksizes = BlockOffsets();
    auto blockContractions = BlockContractions();
#endif

#pragma omp parallel
    {
    auto Cblock = Block(rC,0);

#ifdef ITENSOR_USE_OMP
    int thread_num = omp_get_thread_num();
#endif

    auto Ablocks = std::vector<Block>(A.offsets.size());
    std::transform(A.offsets.begin(),
                   A.offsets.end(),
                   Ablocks.begin(),
                   [](auto pair){return pair.first;});

#pragma omp for schedule(dynamic)
    for(int i=0; i<(int)A.offsets.size(); ++i)
        {
        auto const& Ablock = Ablocks[i];

        //Begin computing elements of Cblock(=destination of this block-block contraction)
        for(auto iA : range(rA))
            if(AtoC[iA] != -1) Cblock[AtoC[iA]] = Ablock[iA];

        //Loop over blocks of B which contract with current block of A
        for(auto const& bio : B.offsets)
            {
            auto do_blocks_contract = true;
            for(auto iA : range(rA))
                {
                auto iB = AtoB[iA];
                if(AtoB[iA] != -1)
                    if(Ablock[iA] != bio.first[iB])
                        {
                        do_blocks_contract = false;
                        break;
                        }
                }
            if(!do_blocks_contract) continue;

            //Finish making Cblock
            for(auto iB : range(rB))
                if(BtoC[iB] != -1) Cblock[BtoC[iB]] = bio.first[iB];

            // Store the current contraction
#ifdef ITENSOR_USE_OMP
            blockContractions_thread[thread_num].push_back(std::make_tuple(Ablock,bio.first,Cblock));
#else
            blockContractions.push_back(std::make_tuple(Ablock,bio.first,Cblock));
#endif

            long blockDim = 1;   //accumulate dim of Indices
            for(auto j : range(order(Cis)))
                {
                auto& J = Cis[j];
                auto i_j = Cblock[j];
                blockDim *= J.blocksize0(i_j);
                }
#ifdef ITENSOR_USE_OMP
            Cblocksizes_thread[thread_num][Cblock] = blockDim;
#else
            Cblocksizes[Cblock] = blockDim;
#endif
            } //for B.offsets
        } //for A.offsets
    }  // omp parallel

#ifdef ITENSOR_USE_OMP
    // Combine blocks from threads
    auto Cblocksizes = BlockOffsets();
    auto blockContractions = BlockContractions();
    for(int thread_num=0; thread_num<num_threads; ++thread_num)
        {
        Cblocksizes.insert(Cblocksizes_thread[thread_num].begin(),
                           Cblocksizes_thread[thread_num].end());
        blockContractions.insert(blockContractions.end(), 
                                 blockContractions_thread[thread_num].begin(),
                                 blockContractions_thread[thread_num].end());
        }
#endif

    auto current_offset = 0;
    for(auto i : Cblocksizes)
        {
        auto current_block = i.first;
        auto current_size = i.second;
        Cblocksizes[current_block] = current_offset;
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
_loopContractedBlocks(QDense<TA> const& A,
                      IndexSet const& Ais,
                      Labels const& Alabel,
                      QDense<TB> const& B,
                      IndexSet const& Bis,
                      Labels const& Blabel,
                      QDense<TC> & C,
                      IndexSet const& Cis,
                      Labels const& Clabel,
                      std::vector<std::tuple<Block,Block,Block>> const& blockContractions,
                      Callable & callback)
    {
    for(auto const& [Ablockind,Bblockind,Cblockind] : blockContractions)
        {
        auto ablock = getBlock(A,Ais,Ablockind);
        auto bblock = getBlock(B,Bis,Bblockind);
        auto cblock = getBlock(C,Cis,Cblockind);
        callback(ablock,Ablockind,
                 bblock,Bblockind,
                 cblock,Cblockind);
        }
    }

long long int inline
contractionCost(Range const& Arange,
                Range const& Brange,
                Labels const& Blabel)
    {
    long long int cost = 1;
    for(decltype(Arange.order()) i = 0; i < Arange.order(); i++)
        cost *= Arange.extent(i);
    for(decltype(Brange.order()) i = 0; i < Brange.order(); i++)
        if(Blabel[i] > 0) cost *= Brange.extent(i);
    return cost;
    }
               

template<typename TA,
         typename TB,
         typename TC,
         typename Callable>
void
_loopContractedBlocksOMP(QDense<TA> const& A,
                         IndexSet const& Ais,
                         Labels const& Alabel,
                         QDense<TB> const& B,
                         IndexSet const& Bis,
                         Labels const& Blabel,
                         QDense<TC> & C,
                         IndexSet const& Cis,
                         Labels const& Clabel,
                         std::vector<std::tuple<Block,Block,Block>> const& blockContractions,
                         Callable & callback)
    {
    auto sortBlockContractions = [](std::tuple<Block,Block,Block> const& t1,
                                    std::tuple<Block,Block,Block> const& t2)
      { 
      return std::get<2>(t1) < std::get<2>(t2);
      };
    auto blockContractionsSorted = blockContractions;
    std::sort(std::begin(blockContractionsSorted),std::end(blockContractionsSorted),
              sortBlockContractions);

    auto nnzblocksC = C.offsets.size();
    auto offset = std::vector<int>(nnzblocksC,0);
    auto nrepeat = std::vector<int>(nnzblocksC,1);

    int nblockC = 0;
    auto ncontractions = blockContractionsSorted.size();

    for(decltype(ncontractions) i = 1; i < ncontractions; i++)
        {
        if(std::get<2>(blockContractionsSorted[i]) == 
           std::get<2>(blockContractionsSorted[i-1]))
            {
            nrepeat[nblockC] += 1;
            }
        else
            {
            nblockC += 1;
            offset[nblockC] = i;
            }
        }

    // 
    // Sort the contraction list by the costs
    // (it doesn't seem to help)
    //

    auto perm = std::vector<int>(offset.size());
    std::iota(perm.begin(), perm.end(), 0);

    /*
    auto costs = std::vector<long long int>(ncontractions);

    for(decltype(ncontractions) i = 0; i < ncontractions; i++)
        {
        auto const& [Ablockind,Bblockind,Cblockind] = blockContractionsSorted[i];
        Range Arange,
              Brange,
              Crange;
        //Construct range objects for aref,bref,cref
        //using IndexDim helper objects
        Arange.init(make_indexdim(Ais,Ablockind));
        Brange.init(make_indexdim(Bis,Bblockind));
        Crange.init(make_indexdim(Cis,Cblockind));
        auto cost = contractionCost(Arange,
                                    Brange,Blabel);
        costs[i] = cost;
        }
      
    auto costs_per_output_block = std::vector<long long int>(nnzblocksC);

    decltype(ncontractions) n = 0;
    for(decltype(nnzblocksC) i = 0; i < nnzblocksC; i++)
      {
      decltype(ncontractions) last_offset_i = offset[i]+nrepeat[i];
      auto& current_cost = costs_per_output_block[i];
      for(decltype(ncontractions) j = offset[i]; j < last_offset_i; j++)
        {
        if(j != n) Error("Wrong contraction plan in QDense contraction");
        current_cost += costs[n];
        n++;
        }
      }

    // Sort offset and nrepeat according to costs_per_output_block
    auto sortby = [&](int i1, int i2) -> bool 
                { 
                return costs_per_output_block[i1] > costs_per_output_block[i2];
                };
    std::sort(perm.begin(), perm.end(), sortby);
    */

#pragma omp parallel for schedule(dynamic)
    for(decltype(nnzblocksC) i = 0; i < nnzblocksC; i++)
      {
      // Contractions that have the same output block
      // location in C are put in the same thread to
      // avoid race conditions
      auto p = perm[i];
      for(auto j = offset[p]; j < offset[p]+nrepeat[p]; j++)
        {
        auto const& [Ablockind,Bblockind,Cblockind] = blockContractionsSorted[j];
        auto ablock = getBlock(A,Ais,Ablockind);
        auto bblock = getBlock(B,Bis,Bblockind);
        auto cblock = getBlock(C,Cis,Cblockind);
        callback(ablock,Ablockind,
                 bblock,Bblockind,
                 cblock,Cblockind);
        }
      }
    }

template<typename TA,
         typename TB,
         typename TC,
         typename Callable>
void
_loopContractedBlocksOMPSingleOutputBlock(QDense<TA> const& A,
                                          IndexSet const& Ais,
                                          Labels const& Alabel,
                                          QDense<TB> const& B,
                                          IndexSet const& Bis,
                                          Labels const& Blabel,
                                          QDense<TC> & C,
                                          IndexSet const& Cis,
                                          Labels const& Clabel,
                                          std::vector<std::tuple<Block,Block,Block>> const& blockContractions,
                                          Callable & callback)
    {
    auto ncontractions = blockContractions.size();

    // Fill with zeros to overwrite garbage
    std::fill(C.begin(),C.end(),0.);

    auto cblock_tot = getBlock(C,Cis,C.offsets.begin()->first);
    auto cblock_size = cblock_tot.size();
    auto cblocks = std::vector<std::vector<TC>>(ncontractions);

#pragma omp parallel for schedule(dynamic)
    for(decltype(ncontractions) i = 0; i < ncontractions; i++)
        {
        auto const& [Ablockind,Bblockind,Cblockind] = blockContractions[i];
        auto ablock = getBlock(A,Ais,Ablockind);
        auto bblock = getBlock(B,Bis,Bblockind);
        auto cblockdata = std::vector<TC>(cblock_size);
        auto cblock = makeDataRange<TC>(cblockdata.data(),cblockdata.size());
        callback(ablock,Ablockind,
                 bblock,Bblockind,
                 cblock,Cblockind);
        cblocks[i] = std::vector<TC>(cblock.data(),cblock.data()+cblock.size());
        }

    for(auto const& cblockdata : cblocks)
        {
        auto cblock = makeDataRange<TC>(cblockdata.data(),cblockdata.size());
        makeRef(cblock_tot,VecRange(cblockdata.size())) += makeRef(cblock,VecRange(cblockdata.size()));
        }
    }

template<typename TA,
         typename TB,
         typename TC,
         typename Callable>
void
loopContractedBlocks(QDense<TA> const& A,
                     IndexSet const& Ais,
                     Labels const& Alabel,
                     QDense<TB> const& B,
                     IndexSet const& Bis,
                     Labels const& Blabel,
                     QDense<TC> & C,
                     IndexSet const& Cis,
                     Labels const& Clabel,
                     std::vector<std::tuple<Block,Block,Block>> const& blockContractions,
                     Callable & callback)
    {
#ifdef ITENSOR_USE_OMP
    if(C.offsets.size() == 1)
        _loopContractedBlocksOMPSingleOutputBlock(A,Ais,Alabel,
                                                  B,Bis,Blabel,
                                                  C,Cis,Clabel,
                                                  blockContractions,
                                                  callback);
    else
        _loopContractedBlocksOMP(A,Ais,Alabel,
                                 B,Bis,Blabel,
                                 C,Cis,Clabel,
                                 blockContractions,
                                 callback);
#else
    _loopContractedBlocks(A,Ais,Alabel,
                          B,Bis,Blabel,
                          C,Cis,Clabel,
                          blockContractions,
                          callback);
#endif
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
            auto ival = aio.first[iA];
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

            auto ablock = makeDataRange(A.data(),aio.second,A.size());

            callback(ablock,aio.first,
                     bblock,Bblockind,
                     cblock,Cblockind);
            } //for couB
        } //for A.offsets
    }


} //namespace itensor

#endif
