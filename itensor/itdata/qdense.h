//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QDENSE_H
#define __ITENSOR_QDENSE_H

#include <vector>
#include "itensor/itdata/task_types.h"
#include "itensor/iqindex.h"
#include "itensor/itdata/itdata.h"
#include "itensor/tensor/types.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/call_rewrite.h"

namespace itensor {

template<typename T>
class QDense;

using QDenseReal = QDense<Real>;
using QDenseCplx = QDense<Cplx>;

struct BlOf
    {
    long block;
    long offset;
    };

template<typename T>
class QDense
    {
    static_assert(not std::is_const<T>::value,
                  "Template argument of QDense must be non-const");
    public:
    using value_type = T;
    using storage_type = std::vector<value_type>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;


    //////////////
    std::vector<BlOf> offsets;
        //^ Block index / data offset pairs.
        //  Assumed that block indices are
        //  in increasing order.

    storage_type store;
        //^ tensor data stored contiguously
    //////////////

    QDense() { }

    QDense(IQIndexSet const& is, 
           QN const& div_);

    //template<typename InputIter>
    //QDense(std::vector<BlOf> const& off,
    //       InputIter && b, InputIter && c)
    //     : offsets(off),
    //       store(b,c)
    //       { }

    template<typename... StoreArgs>
    QDense(std::vector<BlOf> const& off,
           StoreArgs&&... sargs)
         : offsets(off),
           store(std::forward<StoreArgs>(sargs)...)
           { }

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
    getElt(IndexSetT<IQIndex> const& is,
           Indexable const& ind) const;

    template<typename Indexable>
    value_type *
    getElt(IndexSetT<IQIndex> const& is,
           Indexable const& ind)
        {
        const auto& cthis = *this;
        return const_cast<value_type*>(cthis.getElt(is,ind));
        }

    long
    updateOffsets(IndexSetT<IQIndex> const& is,
                  QN const& div);

    };

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
doTask(GetElt<IQIndex>& G, QDenseReal const& d);
Cplx
doTask(GetElt<IQIndex>& G, QDenseCplx const& d);

template<typename T>
void
doTask(SetElt<Real,IQIndex>& S, QDense<T>& d);

void
doTask(SetElt<Cplx,IQIndex>& S, QDenseReal const& d, ManageStore & m);

void
doTask(SetElt<Cplx,IQIndex>& S, QDenseCplx & d);

template<typename T>
Cplx
doTask(SumEls<IQIndex>, QDense<T> const& d);

template<typename T>
void
doTask(Mult<Real> const& M, QDense<T>& d);

void
doTask(Mult<Cplx> const& M, QDense<Real> const& d, ManageStore & m);

void
doTask(Mult<Cplx> const& M, QDense<Cplx> & d);

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
void
doTask(PrintIT<IQIndex>& P, QDense<T> const& d);

auto inline constexpr
doTask(StorageType const& S, QDenseReal const& d) ->StorageType::Type { return StorageType::QDenseReal; }

auto inline constexpr
doTask(StorageType const& S, QDenseCplx const& d) ->StorageType::Type { return StorageType::QDenseCplx; }

template<typename TA, typename TB>
void
doTask(PlusEQ<IQIndex> const& P,
       QDense<TA>      const& A,
       QDense<TB>      const& B,
       ManageStore          & m);

template<typename VA, typename VB>
void
doTask(Contract<IQIndex>& Con,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m);

template<typename VA, typename VB>
void
doTask(NCProd<IQIndex>& P,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m);



// Does a binary search over offsets to see 
// if they contain "blockind"
// If so, return the corresponding data offset,
// otherwise return -1
long
offsetOf(std::vector<BlOf> const& offsets,
         long blockind);

template<typename T>
template<typename Indexable>
T const* QDense<T>::
getElt(IQIndexSet const& is,
       Indexable const& ind) const
    {
    auto r = long(ind.size());
    if(r == 0) return store.data();
#ifdef DEBUG
    if(is.r() != r) 
        {
        printfln("is.r() = %d, ind.size() = %d",is.r(),ind.size());
        Error("Mismatched size of IQIndexSet and elt_ind in get_block");
        }
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
    auto boff = offsetOf(offsets,bind);
    if(boff >= 0)
        {
#ifdef DEBUG
        if(size_t(boff+eoff) >= store.size()) Error("get_elt out of range");
#endif
        return store.data()+boff+eoff;
        }
    return nullptr;
    }

template<typename T>
void
doTask(Order<IQIndex> const& P,
       QDense<T>           & dA);

template<typename T>
void
permuteQDense(Permutation const& P,
             QDense<T>    const& dA,
             IQIndexSet   const& Ais,
             QDense<T>         & dB,
             IQIndexSet   const& Bis);

//template<typename BlockSparseStore, typename Indexable>
//auto
//getBlock(BlockSparseStore & d,
//         IQIndexSet const& is,
//         Indexable const& block_ind)
//    -> decltype(getBlock(std::declval<const BlockSparseStore>(),is,block_ind).cast_away_const())
//    {
//    auto const& cd = d;
//    //ugly but safe, efficient, and avoids code duplication (Meyers, Effective C++)
//    return getBlock(cd,is,block_ind).cast_away_const();
//    }

QN
calcDiv(IQIndexSet const& is, 
        Labels const& block_ind);

//code for doTask(ToITensor...) is in iqtensor.cc
template<typename V>
ITensor
doTask(ToITensor & T, QDense<V> const& d);

template<typename V>
bool
doTask(IsEmpty, QDense<V> const& d) { return d.offsets.empty(); }

} //namespace itensor

#endif

