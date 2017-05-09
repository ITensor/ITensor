//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QCOMBINER_H
#define __ITENSOR_QCOMBINER_H

#include "itensor/itdata/qdense.h"
#include <tuple>

namespace itensor {

class QN;

class QCombiner
    {
    struct BlockRange
        {
        size_t block  = 0,
               start  = 0,
               extent = 0;

        void 
        write(std::ostream& s) const
            { 
            s.write((char*)&block, sizeof(size_t));
            s.write((char*)&start, sizeof(size_t));
            s.write((char*)&extent, sizeof(size_t));
            }
    
        void 
        read(std::istream& s)
            { 
            s.read((char*)&block, sizeof(size_t));
            s.read((char*)&start, sizeof(size_t));
            s.read((char*)&extent, sizeof(size_t));
            }
        };
    public:
    using storage_type = std::vector<BlockRange>;
    using size_type = storage_type::size_type;
    //private:
    Range R_;
    storage_type store_;
    public:

    QCombiner() { }
    
    template<typename IQInds>
    explicit
    QCombiner(IQInds const& cinds)
        { 
        //set up range to sum over all possible
        //blocks that can be formed out of combined inds
        auto RB = RangeBuilder(cinds.size());
        for(auto j : itensor::range(cinds))
            RB.nextIndex(cinds[j].nindex());
        R_ = RB.build();
        store_.resize(area(R_));
        }

    explicit
    QCombiner(Range && range, storage_type && store)
      : R_(std::move(range)),
        store_(std::move(store))
        { }

    Range const&
    range() const { return R_; }

    storage_type const&
    store() const { return store_; }

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

const char*
typeNameOf(QCombiner const& d);

void
read(std::istream& s, QCombiner & dat);

void
write(std::ostream& s, QCombiner const& dat);

Cplx
doTask(GetElt<IQIndex> const& g, QCombiner const& c);

Real inline
doTask(NormNoScale, QCombiner const& d) { return 0; }

void inline
doTask(Conj,QCombiner const& d) { }

template<typename T>
void
doTask(Contract<IQIndex> & C,
       QDense<T>    const& d,
       QCombiner  const& cmb,
       ManageStore       & m);

template<typename T>
void
doTask(Contract<IQIndex> & C,
       QCombiner  const& cmb,
       QDense<T>    const& d,
       ManageStore       & m);

void inline
doTask(PrintIT<IQIndex> & P, 
       QCombiner const& d) { P.s << "QCombiner "; }

auto inline
doTask(StorageType const& S, QCombiner const& d) ->StorageType::Type { return StorageType::QCombiner; }

QN inline
doTask(CalcDiv const& C, QCombiner const& d) { return QN{}; }

bool inline
doTask(CheckComplex, QCombiner const& d) { return false; }

} //namespace itensor

#endif

