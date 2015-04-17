//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDATA_H
#define __ITENSOR_IQTDATA_H

#include "itdata.h"
#include "../iqindex.h"
#include "../indexset.h"
#include "../detail/gcounter.h"

namespace itensor {

template<typename T>
class IQTData : public RegisterData<IQTData<T>>
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

    std::vector<T> data;
        //^ tensor data stored contiguously
    //////////////

    IQTData(const IQIndexSet& is, 
            const QN& Q);

    template<typename Indexable>
    const T*
    getBlock(const IQIndexSet& is,
             const Indexable& block_ind) const;

    template<typename Indexable>
    T*
    getBlock(const IQIndexSet& is,
             const Indexable& block_ind)
        {
        //ugly but safe, efficient, and avoids code duplication (Meyers, Effective C++)
        return const_cast<T*>(static_cast<const IQTData&>(*this).getBlock(is,block_ind));
        }

    template<typename Indexable>
    const T*
    getElt(const IQIndexSet& is,
           const Indexable& ind) const;

    template<typename Indexable>
    T*
    getElt(const IQIndexSet& is,
           const Indexable& ind)
        {
        return const_cast<T*>(static_cast<const IQTData&>(*this).getElt(is,ind));
        }

    long
    updateOffsets(const IQIndexSet& is,
                  const QN& Q);

    virtual
    ~IQTData() { }

    private:


    };

namespace detail {

//function object for calling binaryFind
//on offset vectors below
template<typename T>
struct compBlock
    {
    using BlockOffset = typename IQTData<T>::BlockOffset;
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

}; //namespace detail

template<typename T>
long IQTData<T>::
updateOffsets(const IQIndexSet& is,
              const QN& Q)
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
            const auto& J = is[j];
            auto i = C.i.fast(j);
            blockqn += J.qn(1+i)*J.dir();
            }
        if(blockqn == Q)
            {
            long indstr = 1, //accumulate Index strides
                 ind = 0,
                 totm = 1;   //accumulate area of Indices
            for(int j = 0; j < is.r(); ++j)
                {
                const auto& J = is[j];
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

template<typename T>
IQTData<T>::
IQTData(const IQIndexSet& is, 
        const QN& Q)
    {
    auto totalsize = updateOffsets(is,Q);
    data.assign(totalsize,0);
    }

template<typename T>
template<typename Indexable>
const T* IQTData<T>::
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
    auto res = detail::binaryFind(offsets,ii,detail::compBlock<T>());
    if(res)
        {
        return data.data()+res->offset;
        }
    return nullptr;
    }

template<typename T>
template<typename Indexable>
const T* IQTData<T>::
getElt(const IQIndexSet& is,
       const Indexable& ind) const
    {
    auto r = long(ind.size());
    if(r == 0) return data.data();
#ifdef DEBUG
    if(is.r() != r) Error("Mismatched size of IQIndexSet and elt_ind in get_block");
#endif
    long boff = 0, //block offset
         bstr = 1, //block stride so far
         eoff = 0, //element offset within block
         estr = 1; //element stride
    for(auto i = 0; i < r; ++i)
        {
        const auto& I = is[i];
        long block_ind = 0,
             elt_ind = ind[i];
        while(elt_ind >= I[block_ind].m()) //elt_ind 0-indexed
            {
            elt_ind -= I[block_ind].m();
            ++block_ind;
            }
        boff += block_ind*bstr;
        bstr *= I.nindex();
        eoff += elt_ind*estr;
        estr *= I[block_ind].m();
        }
    //Do a binary search (equal_range) to see
    //if there is a block with block offset "boff"
    auto res = detail::binaryFind(offsets,boff,detail::compBlock<T>());
    if(res)
        {
#ifdef DEBUG
        if(res->offset+eoff >= data.size()) Error("get_elt out of range");
#endif
        return data.data()+res->offset+eoff;
        }
    return nullptr;
    }

}; //namespace itensor

#endif

