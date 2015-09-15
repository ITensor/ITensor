//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDATA_H
#define __ITENSOR_IQTDATA_H

#include <vector>
#include "itensor/itdata/task_types.h"
#include "itensor/iqindex.h"
#include "itensor/itdata/itdata.h"
#include "itensor/tensor/types.h"
#include "itensor/detail/gcounter.h"

namespace itensor {

class ManageStore;
class IQTReal;


QN
calcDiv(IQIndexSet const& is, Label const& block_ind);

// For a block index (0,1,...,Nblocks-1),
// as in the offsets array of an IQTReal,
// computes the zero-indexed location of
// the given block (e.g. 5 -> {1,0,2})
// storing the resulting indices in the 
// container "ind". The jth index of ind
// can range from 0 to is[j].nindex(), 
// such that these indices correspond to
// viewing the IQTReal storage as a 
// "tensor of tensors"
template<typename Container>
void
computeBlockInd(long block,
                IQIndexSet const& is,
                Container& ind)
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

class IQTReal
    {
    public:

    struct BlOf
        {
        long block;
        long offset;
        };

    //////////////
    std::vector<BlOf> offsets;
        //^ Block index / data offset pairs.
        //  Assumed that block indices are
        //  in increasing order.

    std::vector<Real> store;
        //^ tensor data stored contiguously
    //////////////

    IQTReal() { }

    IQTReal(IQIndexSet const& is, 
            QN const& div_);

    explicit operator bool() const { return !store.empty(); }

    Real*
    data() { return store.data(); }

    const Real*
    data() const { return store.data(); }

    size_t
    size() const { return store.size(); }

    template<typename Indexable>
    const Real*
    getElt(IndexSetT<IQIndex> const& is,
           Indexable const& ind) const;

    template<typename Indexable>
    Real*
    getElt(IndexSetT<IQIndex> const& is,
           Indexable const& ind)
        {
        const auto& cthis = *this;
        return const_cast<Real*>(cthis.getElt(is,ind));
        }

    long
    updateOffsets(IndexSetT<IQIndex> const& is,
                  QN const& div);

    };


void inline
write(std::ostream & s, IQTReal const& dat)
    {
    itensor::write(s,dat.offsets);
    itensor::write(s,dat.store);
    }

void inline
read(std::istream & s, IQTReal & dat)
    {
    itensor::read(s,dat.offsets);
    itensor::read(s,dat.store);
    }

void inline
swap(IQTReal & d1,
     IQTReal & d2)
    {
    d1.offsets.swap(d2.offsets);
    d1.store.swap(d2.store);
    }

QN
doTask(CalcDiv const& C, IQTReal const& D);

template <typename F>
void
doTask(ApplyIT<F>& A, IQTReal& d)
    {
    for(auto& elt : d.store)
        elt = A.f(elt);
    }

template <typename F>
void
doTask(VisitIT<F>& V, const IQTReal& d)
    {
    for(const auto& elt : d.store)
        V.f(elt*V.scale_fac);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, IQTReal& d)
    {
    std::generate(d.store.begin(),d.store.end(),G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, const IQTReal& cd, ManageStore& mp)
    {
    Error("Complex version of IQTensor generate not yet supported");
    }


Cplx
doTask(GetElt<IQIndex>& G, const IQTReal& d);

void
doTask(SetElt<Real,IQIndex>& S, IQTReal& d);

//void
//doTask(SetElt<Cplx,IQIndex>& S, IQTReal& d);

void
doTask(MultReal const& M, IQTReal& d);

void
doTask(const PlusEQ<IQIndex>& P,
       IQTReal& A,
       const IQTReal& B);

void
doTask(Contract<IQIndex>& Con,
       const IQTReal& A,
       const IQTReal& B,
       ManageStore& mp);

void
doTask(Conj, const IQTReal& d);

bool inline
doTask(CheckComplex,const IQTReal& d) { return false; }

Real
doTask(NormNoScale, const IQTReal& d);

void
doTask(PrintIT<IQIndex>& P, const IQTReal& d);

void
doTask(Write& W, const IQTReal& d);

IQTReal::BlOf inline
make_blof(long b, long o)
    {
    IQTReal::BlOf B;
    B.block = b;
    B.offset = o;
    return B;
    }

// Does a binary search over offsets to see 
// if they contain "blockind"
// If so, return the corresponding data offset,
// otherwise return -1
long
offsetOf(std::vector<IQTReal::BlOf> const& offsets,
         long blockind);


//
// Helper object for treating
// IQTReal storage as a "tensor of tensors"
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
    { 
    return IndexDim<Indexable>(is,ind); 
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
    auto rA = Ais.r(),
         rB = Bis.r(),
         rC = Cis.r();

    Label AtoB(rA,-1),
          AtoC(rA,-1),
          BtoC(rB,-1);
    for(decltype(rC) ic = 0; ic < rC; ++ic)
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
    for(int ia = 0; ia < rA; ++ia)
    for(int ib = 0; ib < rB; ++ib)
        if(Ais[ia] == Bis[ib])
            {
            AtoB[ia] = ib;
            break;
            }

    detail::GCounter couB(rB);
    Label Ablockind(rA,0),
          Cblockind(rC,0);
    //Loop over blocks of A (labeled by elements of A.offsets)
    for(auto& aio : A.offsets)
        {
        //Reconstruct indices labeling this block of A, put into Ablock
        computeBlockInd(aio.block,Ais,Ablockind);
        //Reset couB to run over indices of B (at first)
        couB.reset();
        for(decltype(rB) ib = 0; ib < rB; ++ib)
            couB.setRange(ib,0,Bis[ib].nindex()-1);
        for(decltype(rA) iA = 0; iA < rA; ++iA)
            {
            auto ival = Ablockind[iA];
            //Restrict couB to be fixed for indices of B contracted with A
            if(AtoB[iA] != -1) couB.setRange(AtoB[iA],ival,ival);
            //Begin computing elements of Cblock(=destination of this block-block contraction)
            if(AtoC[iA] != -1) Cblockind[AtoC[iA]] = ival;
            }
        //Loop over blocks of B which contract with current block of A
        for(;couB.notDone(); ++couB)
            {
            //Check whether B contains non-zero block for this setting of couB
            //TODO: check whether block is present by computing its QN flux,
            //      could be faster than calling getBlock
            auto bblock = getBlock(B,Bis,couB.i);
            if(!bblock) continue;

            //Finish making Cblockind and Bblockind
            Label Bblockind(rB,0);
            for(decltype(rB) ib = 0; ib < rB; ++ib)
                {
                if(BtoC[ib] != -1) Cblockind[BtoC[ib]] = couB.i[ib];
                Bblockind[ib] = couB.i[ib];
                }

            auto cblock = getBlock(C,Cis,Cblockind);
            assert(cblock);

            auto ablock = cData(A.data(),aio.offset,A.size());
            assert(ablock);

            callback(ablock,Ablockind,
                     bblock,Bblockind,
                     cblock,Cblockind);

            } //for couB
        } //for A.offsets
    }

template<typename BlockSparseStore, typename Indexable>
cData
getBlock(BlockSparseStore const& d,
         IQIndexSet const& is,
         Indexable const& block_ind)
    {
    auto r = long(block_ind.size());
    if(r == 0) return cData(d.data(),d.size());
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
    if(boff >= 0) return cData(d.data(),boff,d.size());
    return cData{};
    }

template<typename BlockSparseStore, typename Indexable>
Data
getBlock(BlockSparseStore & d,
         IQIndexSet const& is,
         Indexable const& block_ind)
    {
    auto const& cd = d;
    //ugly but safe, efficient, and avoids code duplication (Meyers, Effective C++)
    return getBlock(cd,is,block_ind).cast_away_const();
    }


template<typename Indexable>
const Real* IQTReal::
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

} //namespace itensor

#endif

