//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTCOMBINER_H
#define __ITENSOR_IQTCOMBINER_H

#include "itensor/itdata/iqtreal.h"
#include <tuple>

namespace itensor {

class QN;

class IQTCombiner
    {
    struct BlockRange
        {
        size_t block  = 0,
               start  = 0,
               extent = 0;
        };
    public:
    using storage_type = std::vector<BlockRange>;
    using size_type = storage_type::size_type;
    private:
    Range R_;
    storage_type store_;
    public:
    
    template<typename IQInds>
    explicit
    IQTCombiner(IQInds const& inds)
        { 
        //set up range to sum over all possible
        //blocks that can be formed out of inds
        auto RB = RangeBuilder(inds.size());
        for(decltype(inds.size()) j = 0; j < inds.size(); ++j)
            RB.nextIndex(inds[j].nindex());
        R_ = RB.build();
        store_.resize(area(R_));
        }

    Range const&
    range() const { return R_; }

    void
    setBlockRange(Range::iterator const& bit,
                  size_type block,
                  size_type start,
                  size_type extent)
        {
        auto& bs = store_.at(bit.offset());
        bs.block = block;
        bs.start = start;
        bs.extent = extent;
        }
                  
    template<typename BlockInd>
    std::tuple<size_type,size_type,size_type>
    getBlockRange(BlockInd const& bind) const
        {
        auto& bs = store_.at(offset(R_,bind));
        return std::make_tuple(bs.block,bs.start,bs.start+bs.extent);
        }

    };

void inline
read(std::istream& s, IQTCombiner & dat) { }

void inline
write(std::ostream& s, IQTCombiner const& dat) { }

Cplx
doTask(GetElt<IQIndex> const& g, IQTCombiner const& c);

Real inline
doTask(NormNoScale, IQTCombiner const& d) { return 0; }

void inline
doTask(Conj,IQTCombiner const& d) { }

void
doTask(Contract<IQIndex> & C,
       IQTReal      const& d,
       IQTCombiner  const& cmb,
       ManageStore       & m);

void
doTask(Contract<IQIndex> & C,
       IQTCombiner  const& cmb,
       IQTReal      const& d,
       ManageStore       & m);

void inline
doTask(PrintIT<IQIndex> & P, 
       IQTCombiner const& d) { P.s << "IQTCombiner "; }

void 
doTask(Write& W, IQTCombiner const& d);

QN inline
doTask(CalcDiv const& C, IQTCombiner const& d) { return QN{}; }

bool inline
doTask(CheckComplex, IQTCombiner const& d) { return false; }

} //namespace itensor

#endif

