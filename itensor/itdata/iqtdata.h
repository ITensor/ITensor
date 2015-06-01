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

class ManagePtr;
class ITCombiner;
class IQTData;


QN
calcDiv(const IQIndexSet& is, const IQTData& D);

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

class IQTData : public RegisterData<IQTData>
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
write(std::ostream& s, const IQTData& dat)
    {
    itensor::write(s,dat.offsets);
    itensor::write(s,dat.data);
    }

void inline
read(std::istream& s, IQTData& dat)
    {
    itensor::read(s,dat.offsets);
    itensor::read(s,dat.data);
    }

//
// Helper object for treating
// IQTData storage as a "tensor of tensors"
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
doTask(ApplyIT<F>& A, IQTData& d)
    {
    for(auto& elt : d.data)
        elt = A.f(elt);
    }

template <typename F>
void
doTask(VisitIT<F>& V, const IQTData& d)
    {
    for(const auto& elt : d.data)
        V.f(elt*V.scale_fac);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, IQTData& d)
    {
    std::generate(d.data.begin(),d.data.end(),G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, const IQTData& cd, ManagePtr& mp)
    {
    Error("Complex version of IQTensor generate not yet supported");
    }


Cplx
doTask(GetElt<IQIndex>& G, const IQTData& d);

void
doTask(SetElt<Real,IQIndex>& S, IQTData& d);

//void
//doTask(SetElt<Cplx,IQIndex>& S, IQTData& d);

void
doTask(MultReal& M, IQTData& d);

void
doTask(const PlusEQ<IQIndex>& P,
       IQTData& A,
       const IQTData& B);

void
doTask(Contract<IQIndex>& Con,
       const IQTData& A,
       const IQTData& B,
       ManagePtr& mp);

void
doTask(Contract<IQIndex>& C,
       const IQTData& d,
       const ITCombiner& cmb,
       ManagePtr& mp);

void
doTask(Contract<IQIndex>& C,
       const ITCombiner& cmb,
       const IQTData& d,
       ManagePtr& mp);

void
doTask(Conj, const IQTData& d);

bool inline
doTask(CheckComplex,const IQTData& d) { return false; }

Real
doTask(const NormNoScale<IQIndex>& N, const IQTData& d);

void
doTask(PrintIT<IQIndex>& P, const IQTData& d);

void
doTask(Write& W, const IQTData& d);

} //namespace itensor

#endif

