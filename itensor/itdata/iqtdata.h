//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDATA_H
#define __ITENSOR_IQTDATA_H

#include <vector>
#include "itensor/matrix/types.h"
#include "itensor/indexset.h"

namespace itensor {

class IQIndex;
class QN;
class IQTData;

QN
calcDiv(const IQIndexSet& is, const IQTData& D);

QN
calcDiv(const IQIndexSet& is, const std::vector<long>& block_ind);

void
inverseBlockInd(long block,
                const IQIndexSet& is,
                std::vector<long>& ind);

class IQTData
    {
    public:

    struct BlockOffset
        {
        long block = 0;
        long offset = 0;
        BlockOffset(long b, long o) : block(b), offset(o) { }
        };

    //////////////
    //Data Members:
    std::vector<BlockOffset> offsets;
        //^ Block index / data offset pairs.
        //Assumed that block indices are
        //in increasing order.

    std::vector<Real> data;
        //^ tensor data stored contiguously
    //////////////

    IQTData(const IQIndexSet& is, 
            const QN& div_);

    IQTData(IQTData&& other) :
        offsets(std::move(other.offsets)),
        data(std::move(other.data))
        { }


    explicit operator bool() const { return !data.empty(); }

    template<typename Indexable>
    const Real*
    getBlock(const IQIndexSet& is,
             const Indexable& block_ind) const;

    template<typename Indexable>
    Real*
    getBlock(const IndexSetT<IQIndex>& is,
             const Indexable& block_ind)
        {
        //ugly but safe, efficient, and avoids code duplication (Meyers, Effective C++)
        return const_cast<Real*>(static_cast<const IQTData&>(*this).getBlock(is,block_ind));
        }

    template<typename Indexable>
    const Real*
    getElt(const IndexSetT<IQIndex>& is,
           const Indexable& ind) const;

    template<typename Indexable>
    Real*
    getElt(const IndexSetT<IQIndex>& is,
           const Indexable& ind)
        {
        return const_cast<Real*>(static_cast<const IQTData&>(*this).getElt(is,ind));
        }

    long
    offsetOf(long blkind) const;

    long
    updateOffsets(const IndexSetT<IQIndex>& is,
                  const QN& div);

    virtual
    ~IQTData() { }

    };

} //namespace itensor

#endif

