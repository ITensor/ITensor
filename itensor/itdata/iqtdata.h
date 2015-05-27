//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDATA_H
#define __ITENSOR_IQTDATA_H

#include <vector>
#include "itensor/matrix/types.h"
#include "itensor/indexset.h"
#include "itensor/iqindex.h"

namespace itensor {

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

    IQTData() { }

    IQTData(const IQIndexSet& is, 
            const QN& div_);


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

template<typename Indexable>
const Real* IQTData::
getBlock(const IQIndexSet& is,
         const Indexable& block_ind) const
    {
    auto r = long(block_ind.size());
    if(r == 0) return data.data();
#ifdef DEBUG
    if(is.r() != r) Error("Mismatched size of IQIndexSet and block_ind in get_block");
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
    auto boff = offsetOf(ii);
    if(boff >= 0)
        {
        return data.data()+boff;
        }
    return nullptr;
    }

template<typename Indexable>
const Real* IQTData::
getElt(const IQIndexSet& is,
       const Indexable& ind) const
    {
    auto r = long(ind.size());
    if(r == 0) return data.data();
#ifdef DEBUG
    if(is.r() != r) Error("Mismatched size of IQIndexSet and elt_ind in get_block");
#endif
    long bind = 0, //block index (total)
         bstr = 1, //block stride so far
         eoff = 0, //element offset within block
         estr = 1; //element stride
    for(auto i = 0; i < r; ++i)
        {
        auto& I = is[i];
        long block_subind = 0,
             elt_subind = ind[i];
        while(elt_subind >= I[block_subind].m()) //elt_ind 0-indexed
            {
            elt_subind -= I[block_subind].m();
            ++block_subind;
            }
        bind += block_subind*bstr;
        bstr *= I.nindex();
        eoff += elt_subind*estr;
        estr *= I[block_subind].m();
        }
    //Do a binary search (equal_range) to see
    //if there is a block with block index "bind"
    auto boff = offsetOf(bind);
    if(boff >= 0)
        {
#ifdef DEBUG
        if(size_t(boff+eoff) >= data.size()) Error("get_elt out of range");
#endif
        return data.data()+boff+eoff;
        }
    return nullptr;
    }

} //namespace itensor

#endif

