//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDIAG_H
#define __ITENSOR_IQTDIAG_H

#include "itensor/itdata/iqtreal.h"

namespace itensor {

class ManageStore;
class ITCombiner;
class IQTReal;


class IQTDiag
    {
    public:

    //////////////
    std::vector<IQTReal::BlOf> offsets;
        //^ Block index / data offset pairs.
        //  Assumed that block indices are
        //  in increasing order.

    std::vector<Real> store;
        //^ *diagonal* tensor data stored contiguously
    //////////////

    IQTDiag() { }

    IQTDiag(IQIndexSet const& is, 
            QN const& div_);

    explicit operator bool() const { return !store.empty(); }

    Real*
    data() { return store.data(); }

    const Real*
    data() const { return store.data(); }

    size_t
    size() const { return store.size(); }

    long
    updateOffsets(IQIndexSet const& is,
                  QN const& div);

    };


template<typename Indexable>
const Real*
getElt(IQTDiag const& D,
       IQIndexSet const& is,
       Indexable const& ind)
    {
    auto r = is.r();
#ifdef DEBUG
    if(is.r() != decltype(r)(ind.size())) 
        {
        printfln("is.r() = %d, ind.size() = %d",is.r(),ind.size());
        Error("Wrong number of indices passed to getElt");
        }
#endif
    if(r == 0) return D.data();
    long bind = 0, //block index (total)
         bstr = 1; //block stride so far
    auto last_elt_subind = ind[0];
    for(decltype(r) i = 0; i < r; ++i)
        {
        auto& I = is[i];
        long block_subind = 0,
             elt_subind = ind[i];
        while(elt_subind >= I[block_subind].m()) //elt_subind 0-indexed
            {
            elt_subind -= I[block_subind].m();
            ++block_subind;
            }
        if(i != 0 && elt_subind != last_elt_subind) return nullptr;

        last_elt_subind = elt_subind;
        bind += block_subind*bstr;
        bstr *= I.nindex();
        }
    //Do a binary search (equal_range) to see
    //if there is a block with block index "bind"
    auto boff = offsetOf(D.offsets,bind);
    if(boff != -1)
        {
        auto eoff = last_elt_subind;
#ifdef DEBUG
        if(size_t(boff+eoff) >= D.store.size()) Error("get_elt out of range");
#endif
        return D.data()+boff+eoff;
        }
    return nullptr;
    }

template<typename Indexable>
Real*
getElt(IQTDiag & D,
       IndexSetT<IQIndex> const& is,
       Indexable const& ind)
    {
    return const_cast<Real*>(getElt(D,is,ind));
    }

void inline
write(std::ostream & s, IQTDiag const& dat)
    {
    itensor::write(s,dat.offsets);
    itensor::write(s,dat.store);
    }

void inline
read(std::istream & s, IQTDiag & dat)
    {
    itensor::read(s,dat.offsets);
    itensor::read(s,dat.store);
    }
 
Cplx
doTask(GetElt<IQIndex>& G, IQTDiag const& D);

QN
doTask(CalcDiv const& C, IQTDiag const& d);

template<typename F>
void
doTask(ApplyIT<F> & A, IQTDiag & d)
    {
    for(auto& elt : d.store)
        elt = A.f(elt);
    }

template<typename F>
void
doTask(VisitIT<F> & V, IQTDiag const& d)
    {
    for(auto& elt : d.store)
        V.f(elt*V.scale_fac);
    }

template<typename F>
void
doTask(GenerateIT<F,Real> & G, IQTDiag & d)
    {
    std::generate(d.store.begin(),d.store.end(),G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx> & G, IQTDiag const& cd, ManageStore & mp)
    {
    Error("Complex version of IQTensor generate not yet supported");
    }

void
doTask(MultReal & M, IQTDiag & d);

void
doTask(Contract<IQIndex>& Con,
       IQTReal const& A,
       IQTDiag const& B,
       ManageStore & mp);

void
doTask(Contract<IQIndex>& Con,
       IQTDiag const& A,
       IQTReal const& B,
       ManageStore & mp);

void
doTask(Conj, IQTDiag const& d);

bool inline
doTask(CheckComplex, IQTDiag const& d) { return false; }

Real
doTask(NormNoScale, IQTDiag const& d);

void
doTask(PrintIT<IQIndex> & P, IQTDiag const& d);

void
doTask(Write & W, IQTReal const& d);

} //namespace itensor

#endif

