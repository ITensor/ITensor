//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/iqtdata.h"

#include "itensor/iqindex.h"
#include "itensor/indexset.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/algs.h"
#include "itensor/util/count.h"

using std::vector;

namespace itensor {

//function object for calling binaryFind
//on offset vectors below
struct compBlock
    {
    using BlockOffset = typename IQTData::BlockOffset;
    bool
    operator()(const BlockOffset& bo1,
               const BlockOffset& bo2) const
        { return bo1.block < bo2.block; }
    bool
    operator()(const BlockOffset& bo, long blk) const        
        { return bo.block < blk; }
    bool
    operator()(long blk, const BlockOffset& bo) const 
        { return blk < bo.block; }
    };

QN
calcDiv(const IQIndexSet& is, const vector<long>& block_ind)
    {
    QN div;
    for(auto i : count(is.r())) { div += is[i].dir()*is[i].qn(1+block_ind[i]); }
    return div;
    }

void
inverseBlockInd(long block,
                const IQIndexSet& is,
                vector<long>& ind)
    {
    auto r = int(ind.size());
    assert(r == is.r());
    for(int j = 0; j < r-1; ++j)
        {
        ind[j] = block % is[j].nindex();
        block = (block-ind[j])/is[j].nindex();
        }
    ind[r-1] = block;
    }

QN
calcDiv(const IQIndexSet& is, const IQTData& D)
    {
#ifdef DEBUG
    if(D.offsets.empty()) Error("Default constructed IQTData in calcDiv");
#endif
    auto b = D.offsets.front().block;
    QN div;
    auto r = long(is.r());
    for(long j = 0; j < r-1; ++j)
        {
        auto& J = is[j];
        auto Ij = b % J.nindex();
        div += J.dir()*J.qn(1+Ij);
        b = (b-Ij)/J.nindex();
        }
    div += is[r-1].dir()*is[r-1].qn(1+b);
    return div;
    }

IQTData::
IQTData(const IQIndexSet& is, 
        const QN& div)
    {
    auto totalsize = updateOffsets(is,div);
    data.assign(totalsize,0);
    }

long IQTData::
updateOffsets(const IQIndexSet& is,
              const QN& div)
    {
    offsets.clear();

    if(is.r()==0)
        {
        offsets.emplace_back(0,0);
        return 1;
        }

    detail::GCounter C(0,is.r()-1,0);
    for(int j = 0; j < is.r(); ++j) 
        C.setInd(j,0,is[j].nindex()-1);

    long totalsize = 0;
    for(; C.notDone(); ++C)
        {
        QN blockqn;
        for(int j = 0; j < is.r(); ++j)
            {
            auto& J = is[j];
            blockqn += J.qn(1+C.i[j])*J.dir();
            }
        if(blockqn == div)
            {
            long indstr = 1, //accumulate Index strides
                 ind = 0,
                 totm = 1;   //accumulate area of Indices
            for(int j = 0; j < is.r(); ++j)
                {
                auto& J = is[j];
                auto i_j = C.i[j];
                ind += i_j*indstr;
                indstr *= J.nindex();
                totm *= J[i_j].m();
                }
            offsets.emplace_back(ind,totalsize);
            totalsize += totm;
            }
        }
    return totalsize;
    }




long IQTData::
offsetOf(long blkind) const
    {
    auto blk = detail::binaryFind(offsets,blkind,compBlock());
    if(blk) return blk->offset;
    return -1;
    }

} //namespace itensor

