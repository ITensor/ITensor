//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
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
    IQIndexSet const& is_;
    Indexable const& ind_;
    public:

    IndexDim(IQIndexSet const& is,
             Indexable const& ind)
      : is_(is),
        ind_(ind)
        { }

    size_t
    size() const { return is_.r(); }

    size_t
    operator[](size_t j) const { return (is_[j])[ind_[j]].m(); }
    };

template<typename Indexable>
auto
make_indexdim(IQIndexSet const& is, Indexable const& ind) 
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
// can range from 0 to is[j].nindex(), 
// such that these indices correspond to
// viewing the QDense storage as a 
// "tensor of tensors"
template<typename Container>
void
computeBlockInd(long block,
                IQIndexSet const& is,
                Container & ind)
    {
    using size_type = decltype(ind.size());
    size_type r = ind.size();
    assert(r == size_type(is.r()));
    for(size_type j = 0; j < r-1; ++j)
        {
        auto res = std::ldiv(block,is[j].nindex());
        ind[j] = res.rem;
        block = res.quot;
        }
    ind[r-1] = block;
    }

template<typename BlockSparse, typename Indexable>
auto
getBlock(BlockSparse & d,
         IQIndexSet const& is,
         Indexable const& block_ind)
    -> stdx::if_compiles_return<decltype(makeDataRange(d.data(),d.size())),decltype(d.offsets)>
    {
    auto r = long(block_ind.size());
    if(r == 0) return makeDataRange(d.data(),d.size());
#ifdef DEBUG
    if(is.r() != r) Error("Mismatched size of IQIndexSet and block_ind in getBlock");
#endif
    long ii = 0;
    for(auto i = r-1; i > 0; --i)
        {
        ii += block_ind[i];
        ii *= is[i-1].nindex();
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
                     IQIndexSet const& Ais,
                     BlockSparseB const& B,
                     IQIndexSet const& Bis,
                     BlockSparseC & C,
                     IQIndexSet const& Cis,
                     Callable & callback)
    {
    auto rA = Ais.r();
    auto rB = Bis.r();
    auto rC = Cis.r();

    auto AtoB = IntArray(rA,-1);
    auto AtoC = IntArray(rA,-1);
    auto BtoC = IntArray(rB,-1);
    for(auto ic : range(rC))
        {
        auto j = findindex(Ais,Cis[ic]);
        if(j >= 0)
            {
            AtoC[j] = ic;
            }
        else
            {
            j = findindex(Bis,Cis[ic]);
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
        TIMER_START(19)
        //Reconstruct indices labeling this block of A, put into Ablock
        //TODO: optimize away need to call computeBlockInd by
        //      storing block indices directly in QDense
        //      Taking 10% of running time in S=1 N=100 DMRG tests (maxm=100)
        computeBlockInd(aio.block,Ais,Ablockind);
        //Reset couB to run over indices of B (at first)
        couB.reset();
        for(auto iB : range(rB))
            {
            couB.setRange(iB,0,Bis[iB].nindex()-1);
            }
        for(auto iA : range(rA))
            {
            auto ival = Ablockind[iA];
            //Restrict couB to be fixed for indices of B contracted with A
            if(AtoB[iA] != -1) couB.setRange(AtoB[iA],ival,ival);
            //Begin computing elements of Cblock(=destination of this block-block contraction)
            if(AtoC[iA] != -1) Cblockind[AtoC[iA]] = ival;
            }
        TIMER_STOP(19)
        //Loop over blocks of B which contract with current block of A
        for(;couB.notDone(); ++couB)
            {
            TIMER_START(19)
            //START_TIMER(33)
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
            TIMER_STOP(19)

            callback(ablock,Ablockind,
                     bblock,Bblockind,
                     cblock,Cblockind);

            } //for couB
        } //for A.offsets
    }


} //namespace itensor

#endif
