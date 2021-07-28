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

        std::ostream&
        print(std::ostream & s) const
            {
            s << "Block: " << block << ", Start: " << start << ", Extent: " << extent << "\n";
            return s;
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
    
    template<typename Inds>
    explicit
    QCombiner(Inds const& cinds)
        { 
        //set up range to sum over all possible
        //blocks that can be formed out of combined inds
        auto RB = RangeBuilder(cinds.size());
        for(auto j : itensor::range(cinds))
            RB.nextIndex(cinds[j].nblock());
        R_ = RB.build();
        store_.resize(dim(R_));
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

std::ostream&
operator<<(std::ostream& s, QCombiner const& dat);

Cplx
doTask(GetElt const& g, QCombiner const& c);

Real inline
doTask(NormNoScale, QCombiner const& d) { return 0; }

void inline
doTask(Conj,QCombiner const& d) { }

template<typename T>
void
doTask(Contract & C,
       QDense<T>    const& d,
       QCombiner  const& cmb,
       ManageStore       & m);

template<typename T>
void
doTask(Contract & C,
       QCombiner  const& cmb,
       QDense<T>    const& d,
       ManageStore       & m);

void inline
doTask(PrintIT & P, 
       QCombiner const& d) { P.s << "QCombiner "; }

auto inline
doTask(StorageType const& S, QCombiner const& d) ->StorageType::Type { return StorageType::QCombiner; }

QN inline
doTask(CalcDiv const& C, QCombiner const& d) { return QN{}; }

bool inline
doTask(CheckComplex, QCombiner const& d) { return false; }

} //namespace itensor

#endif

