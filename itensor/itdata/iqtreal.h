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

namespace itensor {

class ManageStore;
class ITCombiner;
class IQTReal;


QN
calcDiv(const IQIndexSet& is, const IQTReal& D);

QN
calcDiv(const IQIndexSet& is, const Label& block_ind);

template<typename Container>
void
inverseBlockInd(long block,
                const IQIndexSet& is,
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
        //BlockOffset(long b, long o) : block(b), offset(o) { }
        //BlockOffset() { }
        };

    //////////////
    //Data Members:
    std::vector<BlOf> offsets;
        //^ Block index / data offset pairs.
        //Assumed that block indices are
        //in increasing order.

    std::vector<Real> data;
        //^ tensor data stored contiguously
    //////////////

    IQTReal() { }

    IQTReal(const IQIndexSet& is, 
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
        return const_cast<Real*>(static_cast<const IQTReal&>(*this).getBlock(is,block_ind));
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
        return const_cast<Real*>(static_cast<const IQTReal&>(*this).getElt(is,ind));
        }

    long
    offsetOf(long blkind) const;

    long
    updateOffsets(const IndexSetT<IQIndex>& is,
                  const QN& div);

    virtual
    ~IQTReal() { }

    };

template<typename Indexable>
const Real* IQTReal::
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
const Real* IQTReal::
getElt(const IQIndexSet& is,
       const Indexable& ind) const
    {
    auto r = long(ind.size());
    if(r == 0) return data.data();
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


void inline
write(std::ostream& s, const IQTReal& dat)
    {
    itensor::write(s,dat.offsets);
    itensor::write(s,dat.data);
    }

void inline
read(std::istream& s, IQTReal& dat)
    {
    itensor::read(s,dat.offsets);
    itensor::read(s,dat.data);
    }

//
// Helper object for treating
// IQTReal storage as a "tensor of tensors"
//
template<typename Indexable>
class IndexDim
    {
    const IQIndexSet& is_;
    const Indexable& ind_;
    public:

    IndexDim(const IQIndexSet& is,
             const Indexable& ind)
      : is_(is),
        ind_(ind)
        { }

    size_t
    size() const { return is_.r(); }

    size_t
    operator[](size_t j) const { return (is_[j])[ind_[j]].m(); }
    };

template<typename Indexable>
IndexDim<Indexable>
make_indexdim(const IQIndexSet& is, const Indexable& ind) 
    { 
    return IndexDim<Indexable>(is,ind); 
    }


template <typename F>
void
doTask(ApplyIT<F>& A, IQTReal& d)
    {
    for(auto& elt : d.data)
        elt = A.f(elt);
    }

template <typename F>
void
doTask(VisitIT<F>& V, const IQTReal& d)
    {
    for(const auto& elt : d.data)
        V.f(elt*V.scale_fac);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, IQTReal& d)
    {
    std::generate(d.data.begin(),d.data.end(),G.f);
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
doTask(MultReal& M, IQTReal& d);

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
doTask(Contract<IQIndex>& C,
       const IQTReal& d,
       const ITCombiner& cmb,
       ManageStore& m);

void
doTask(Contract<IQIndex>& C,
       const ITCombiner& cmb,
       const IQTReal& d,
       ManageStore& m);

void
doTask(Conj, const IQTReal& d);

bool inline
doTask(CheckComplex,const IQTReal& d) { return false; }

Real
doTask(NormNoScale, const IQTReal& d);

void
doTask(PrintIT<IQIndex>& P, const IQTReal& d);

void
doTask(PrintIT<IQIndex>& P, const ITCombiner& d);

void
doTask(Write& W, const IQTReal& d);

} //namespace itensor

#endif

