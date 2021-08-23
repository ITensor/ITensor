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
#ifndef __ITENSOR_QDENSE_H
#define __ITENSOR_QDENSE_H

#include <vector>
#include "itensor/itdata/task_types.h"
#include "itensor/itdata/itdata.h"
#include "itensor/tensor/types.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/call_rewrite.h"

namespace itensor {

template<typename T>
class QDense;

using QDenseReal = QDense<Real>;
using QDenseCplx = QDense<Cplx>;

template<typename T>
class QDense
    {
    static_assert(not std::is_const<T>::value,
                  "Template argument of QDense must be non-const");
    public:
    using value_type = T;
    using storage_type = vector_no_init<value_type>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;


    //////////////
    BlockOffsets offsets;
        //^ Block index / data offset pairs.
        //  Assumed that block indices are
        //  in increasing order.

    storage_type store;
        //^ tensor data stored contiguously
    //////////////

    QDense() { }

    QDense(IndexSet const& is, 
           QN       const& div);

    // Constructor taking a list of block labels
    // instead of QN divergence
    QDense(IndexSet const& is,
           Blocks   const& blocks);

    //template<typename InputIter>
    //QDense(BlockOffsets const& off,
    //       InputIter && b, InputIter && c)
    //     : offsets(off),
    //       store(b,c)
    //       { }

    QDense(BlockOffsets const& off,
           size_t size)
         : offsets(off),
           store(size)
           {
           std::fill(store.begin(),store.end(),0.);
           }

    QDense(BlockOffsets const& off,
           std::vector<value_type> const& v)
         : offsets(off),
           store(v.begin(),v.end())
           {
           }

    template<typename InputIterator>
    QDense(BlockOffsets const& off,
           InputIterator begin, InputIterator end)
     : offsets(off),
       store(begin,end)
       {
       }

    QDense(UndefInitializer,
           std::vector<BlOf> const& off,
           size_t size)
         : offsets(off),
           store(size)
           {
           }

    // TODO: any need for this generic constructor?
    //template<typename... StoreArgs>
    //QDense(BlockOffsets const& off,
    //       StoreArgs&&... sargs)
    //     : offsets(off),
    //       store(std::forward<StoreArgs>(sargs)...)
    //       {
    //       std::fill(store.begin(),store.end(),0.);
    //       }

    explicit operator bool() const { return !store.empty() && !offsets.empty(); }

    value_type *
    data() { return store.data(); }

    value_type const*
    data() const { return store.data(); }

    size_t
    size() const { return store.size(); }

    iterator
    begin() { return store.begin(); }

    iterator
    end() { return store.end(); }

    const_iterator
    begin() const { return store.begin(); }

    const_iterator
    end() const { return store.end(); }

    template<typename Indexable>
    value_type const*
    getElt(IndexSet const& is,
           Indexable const& ind) const
        {
        return std::get<0>(getEltBlockOffset(is,ind));
        }

    template<typename Indexable>
    value_type *
    getElt(IndexSet const& is,
           Indexable const& ind)
        {
        const auto& cthis = *this;
        return const_cast<value_type*>(cthis.getElt(is,ind));
        }

    template<typename Indexable>
    std::tuple<value_type const*,Block,long>
    getEltBlockOffset(IndexSet const& is,
                      Indexable const& ind) const;

    template<typename Indexable>
    std::tuple<value_type*,Block,long>
    getEltBlockOffset(IndexSet const& is,
                      Indexable const& ind)
        {
        const auto& cthis = *this;
        auto eltblockoffset = cthis.getEltBlockOffset(is,ind);
        return std::make_tuple(const_cast<value_type*>(std::get<0>(eltblockoffset)),
                               std::get<1>(eltblockoffset),
                               std::get<2>(eltblockoffset));
        }

    long
    updateOffsets(IndexSet const& is,
                  QN const& div);

    long
    updateOffsets(IndexSet const& is,
                  Blocks   const& blocks);

    // Insert a new block into the QDense storage
    // and fill with the specified value.
    // Return the offset of the new block.
    long
    insertBlock(IndexSet const& is,
                Block    const& block,
                T val = 0);

    };

std::ostream&
operator<<(std::ostream & s, BlOf const& t);

std::ostream&
operator<<(std::ostream & s, BlockOffsets const& offsets);

std::ostream&
operator<<(std::ostream & s, Blocks const& offsets);

template<typename T>
std::ostream&
operator<<(std::ostream & s, QDense<T> const& t);

const char*
typeNameOf(QDenseReal const& d);
const char*
typeNameOf(QDenseCplx const& d);

//QDenseCplx inline
//makeCplx(QDenseReal const& DR)
//    {
//    auto DC = QDenseCplx{};
//    DC.offsets = DR.offsets;
//    }

template<typename T>
bool constexpr
isReal(QDense<T> const& t) { return std::is_same<T,Real>::value; }

template<typename T>
bool constexpr
isCplx(QDense<T> const& t) { return std::is_same<T,Cplx>::value; }

Data inline
realData(QDenseReal & d) { return Data(d.data(),d.size()); }

Datac inline
realData(QDenseReal const& d) { return Datac(d.data(),d.size()); }

Data inline
realData(QDenseCplx & d) { return Data(reinterpret_cast<Real*>(d.data()),2*d.size()); }

Datac inline
realData(QDenseCplx const& d) { return Datac(reinterpret_cast<const Real*>(d.data()),2*d.size()); }

template<typename T>
void
write(std::ostream & s, QDense<T> const& dat)
    {
    itensor::write(s,dat.offsets);
    itensor::write(s,dat.store);
    }

template<typename T>
void
read(std::istream & s, QDense<T> & dat)
    {
    itensor::read(s,dat.offsets);
    itensor::read(s,dat.store);
    }

template<typename T>
void
swap(QDense<T> & d1,
     QDense<T> & d2)
    {
    d1.offsets.swap(d2.offsets);
    d1.store.swap(d2.store);
    }

template<typename T>
QN
doTask(CalcDiv const& C, QDense<T> const& D);

template <typename F, typename T>
void
doTask(ApplyIT<F>& A, QDense<T>& d)
    {
    if(!d) throw ITError("Empty storage in QDense apply");
    //for(auto& elt : d.store)
    //    elt = A.f(elt);
    for(auto& el : d.store) A(el);
    }

template <typename F, typename T>
void
doTask(VisitIT<F>& V, QDense<T> const& d)
    {
    if(!d) throw ITError("Empty storage in QDense visit");
    for(const auto& el : d.store)
        {
        detail::call<void>(V.f,el*V.scale_fac);
        }
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, QDenseReal & D)
    {
    stdx::generate(D,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, QDenseCplx const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<QDenseReal>(D.offsets,D.size());
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, QDenseReal const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<QDenseCplx>(D.offsets,D.size());
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, QDenseCplx & D)
    {
    stdx::generate(D,G.f);
    }


Cplx
doTask(GetElt& G, QDenseReal const& d);
Cplx
doTask(GetElt& G, QDenseCplx const& d);

template<typename T>
void
doTask(SetElt<Real>& S, QDense<T>& d);

void
doTask(SetElt<Cplx>& S, QDenseReal const& d, ManageStore & m);

void
doTask(SetElt<Cplx>& S, QDenseCplx & d);

template<typename T>
Cplx
doTask(SumEls, QDense<T> const& d);

template<typename T>
void
doTask(Mult<Real> const& M, QDense<T>& d);

void
doTask(Mult<Cplx> const& M, QDense<Real> const& d, ManageStore & m);

void
doTask(Mult<Cplx> const& M, QDense<Cplx> & d);

void
doTask(MakeCplx const&, QDense<Cplx> & d);
     
void
doTask(MakeCplx const&, QDense<Real> const& d, ManageStore & m);

template<typename T>
void
doTask(Fill<T> const& F, QDense<T> & d);

template<typename FT, typename DT,
         class=stdx::require_not<std::is_same<FT,DT>> >
void
doTask(Fill<FT> const& F, QDense<DT> const& d, ManageStore & m);


void inline
doTask(Conj, QDenseReal const& d) { }

void
doTask(Conj, QDenseCplx & d);

template<typename T>
bool constexpr
doTask(CheckComplex, QDense<T> const& d) { return isCplx(d); }

void inline
doTask(TakeReal, QDenseReal const& d) { }

void
doTask(TakeReal, QDenseCplx const& d, ManageStore & m);

void
doTask(TakeImag, QDenseReal & d);

void
doTask(TakeImag, QDenseCplx const& d, ManageStore & m);

template<typename T>
Real
doTask(NormNoScale, QDense<T> const& D);

template<typename T>
int
doTask(NNZBlocks, QDense<T> const& D);

template<typename T>
long
doTask(NNZ, QDense<T> const& D);

template<typename T>
void
doTask(PrintIT& P, QDense<T> const& d);

auto inline constexpr
doTask(StorageType const& S, QDenseReal const& d) ->StorageType::Type { return StorageType::QDenseReal; }

auto inline constexpr
doTask(StorageType const& S, QDenseCplx const& d) ->StorageType::Type { return StorageType::QDenseCplx; }

template<typename TA, typename TB>
void
doTask(PlusEQ const& P,
       QDense<TA>      const& A,
       QDense<TB>      const& B,
       ManageStore          & m);

template<typename VA, typename VB>
void
doTask(Contract& Con,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m);

//TODO: complete implementation
//template<typename VA, typename VB>
//void
//doTask(Contract& Con,
//       Dense<VA> const& A,
//       QDense<VB> const& B,
//       ManageStore& m);

template<typename VA, typename VB>
void
doTask(NCProd& P,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m);

// From an indexset and a QN divergence,
// get the list of block-offsets and the
// size of the storage
std::tuple<BlockOffsets,long>
getBlockOffsets(IndexSet const& is,
                QN       const& div);

// Does a binary search over offsets to see 
// if they contain "blockind"
// If so, return the corresponding data offset,
// otherwise return -1
long
offsetOf(BlockOffsets const& offsets,
         Block const& blockind);

int
offsetOfLoc(BlockOffsets const& offsets,
            Block        const& blockind);

template<typename T>
template<typename Indexable>
std::tuple<T const*,Block,long> QDense<T>::
getEltBlockOffset(IndexSet const& is,
                  Indexable const& ind) const
    {
    auto r = long(ind.size());
    if(r == 0) return std::make_tuple(store.data(),Block(0),0);
#ifdef DEBUG
    if(is.order() != r) 
        {
        printfln("is.order() = %d, ind.size() = %d",is.order(),ind.size());
        Error("Mismatched size of IndexSet and elt_ind in get_block");
        }
#endif
    long eoff = 0, //element offset within block
         estr = 1; //element stride
    auto block = Block(r);
    for(auto i = 0; i < r; ++i)
        {
        auto& I = is[i];
        long block_subind = 0,
             elt_subind = ind[i];
        while(elt_subind >= I.blocksize0(block_subind)) //elt_ind 0-indexed
            {
            elt_subind -= I.blocksize0(block_subind);
            ++block_subind;
            }
        block[i] = block_subind;
        eoff += elt_subind*estr;
        estr *= I.blocksize0(block_subind);
        }
    //Do a binary search (equal_range) to see
    //if there is a block with block index "bind"
    auto boff = offsetOf(offsets,block);
    if(boff >= 0)
        {
#ifdef DEBUG
        if(size_t(boff+eoff) >= store.size()) Error("get_elt out of range");
#endif
        return std::make_tuple(store.data()+boff+eoff,block,eoff);
        }
    return std::make_tuple(nullptr,block,eoff);
    }

template<typename T>
long QDense<T>::
insertBlock(IndexSet const& is,
            Block    const& block,
            T val)
    {
    // Get the block dimension
    int blockdim = 1;
    for(auto i : range(is.order()))
        blockdim *= is[i].blocksize0(block[i]);

    // Find where to insert the new block
    // TODO: optimize with a binary search
    size_t insert_loc = 0;
    for(auto const& bof : offsets)
        {
        if(block > bof.block) insert_loc++;
        else break;
        }

    // Get the offset of the new block
    long new_offset;
    if(insert_loc >= offsets.size())
        new_offset = store.size();
    else
        new_offset = offsets[insert_loc].offset;

    // Insert the specified value into the storage
    store.insert(store.begin()+new_offset,blockdim,val);

    // Shift the offsets by the new block dimension
    for(auto i = insert_loc; i < offsets.size(); i++)
        offsets[i].offset += blockdim;

    // Insert the block and offset into the block-offsets list
    offsets.insert(offsets.begin()+insert_loc,make_blof(block,new_offset));
    return new_offset;
    }

template<typename T>
void
doTask(Order const& P,
       QDense<T>           & dA);

template<typename T>
void
permuteQDense(Permutation const& P,
              QDense<T>    const& dA,
              IndexSet   const& Ais,
              QDense<T>         & dB,
              IndexSet   const& Bis);

//template<typename BlockSparseStore, typename Indexable>
//auto
//getBlock(BlockSparseStore & d,
//         IndexSet const& is,
//         Indexable const& block_ind)
//    -> decltype(getBlock(std::declval<const BlockSparseStore>(),is,block_ind).cast_away_const())
//    {
//    auto const& cd = d;
//    //ugly but safe, efficient, and avoids code duplication (Meyers, Effective C++)
//    return getBlock(cd,is,block_ind).cast_away_const();
//    }

QN
calcDiv(IndexSet const& is, 
        Block const& block_ind);

////code for doTask(ToITensor...) is in iqtensor.cc
//template<typename V>
//ITensor
//doTask(ToITensor & T, QDense<V> const& d);

template<typename V>
bool
doTask(IsEmpty, QDense<V> const& d) { return d.offsets.empty(); }

template<typename V>
TenRef<Range,V>
doTask(GetBlock<V> const& G, QDense<V> & d);

template<typename T>
bool
doTask(IsDense,
       QDense<T> const& d);

template<typename V>
void
doTask(RemoveQNs & T, 
       QDense<V> const& d,
       ManageStore & m);

#ifdef ITENSOR_USE_HDF5
template<typename V>
void
h5_write(h5::group parent, std::string const& name, QDense<V> const& D);

template<typename V>
void
h5_read(h5::group parent, std::string const& name, QDense<V> & D);
#endif

} //namespace itensor

#endif

