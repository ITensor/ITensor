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


template<typename T>
Cplx
doTask(GetElt<IQIndex>& G, QDense<T> const& d);

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


BlOf inline
make_blof(long b, long o)
    {
    BlOf B;
    B.block = b;
    B.offset = o;
    return B;
    }

// Does a binary search over offsets to see 
// if they contain "blockind"
// If so, return the corresponding data offset,
// otherwise return -1
long
offsetOf(std::vector<BlOf> const& offsets,
         long blockind);


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

QN
calcDiv(IQIndexSet const& is, 
        Labels const& block_ind);

//code for doTask(ToITensor...) is in iqtensor.cc
template<typename V>
ITensor
doTask(ToITensor & T, QDense<V> const& d);

} //namespace itensor

#endif

