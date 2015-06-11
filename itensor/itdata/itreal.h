//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITREAL_H
#define __ITENSOR_ITREAL_H

#include "itensor/itdata/task_types.h"
#include "itensor/util/readwrite.h"
#include "itensor/detail/call_rewrite.h"
#include "itensor/itdata/itdata.h"

namespace itensor {


class ITReal
    {
    public:
    using storage_type = std::vector<Real>;
    using size_type = storage_type::size_type;
    using iterator = storage_type::iterator;
    using const_iterator = storage_type::const_iterator;
    using value_type = Real;

    //
    // Data members
    //

    storage_type store;

    //
    // Constructors
    //

    ITReal() { }

    ITReal(size_t size) : store(size) { }

    ITReal(size_t size, Real val) : store(size,val) { }

    template<typename InputIterator>
    ITReal(InputIterator b, InputIterator e) : store(b,e) { }

    ITReal(storage_type&& data) : store(std::move(data)) { }

    //
    //std container like methods
    //

    Real&
    operator[](size_type i) 
        { 
#ifdef DEBUG
        return store.at(i);
#else
        return store[i]; 
#endif
        }
    const Real&
    operator[](size_type i) const 
        {
#ifdef DEBUG
        return store.at(i);
#else
        return store[i]; 
#endif
        }

    size_type
    size() const { return store.size(); }
    bool
    empty() const { return store.empty(); }

    Real*
    data() { return store.data(); }
    const Real*
    data() const { return store.data(); }
    
    const_iterator
    cbegin() const { return store.cbegin(); }

    const_iterator
    cend() const { return store.cend(); }

    const_iterator
    begin() const { return store.begin(); }

    const_iterator
    end() const { return store.end(); }

    iterator
    begin() { return store.begin(); }

    iterator
    end() { return store.end(); }
    };

void inline
read(std::istream& s, ITReal& dat)
    {
    itensor::read(s,dat.store);
    }

void inline
write(std::ostream& s, const ITReal& dat)
    {
    itensor::write(s,dat.store);
    }

template<typename F>
void
doTask(ApplyIT<F>& A, ITReal& d)
    { 
    for(auto& elt : d) elt = detail::call<Real>(A.f,elt);
    }

template<typename F>
void
doTask(VisitIT<F>& V, const ITReal& d)
    { 
    for(auto& elt : d) detail::call<void>(V.f,V.scale_fac * elt);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, ITReal& d)
    { 
    std::generate(d.begin(),d.end(),G.f);
    }

Cplx 
doTask(const GetElt<Index>& g, const ITReal& d);

void
doTask(const SetElt<Real,Index>& s, ITReal& d);

void
doTask(const SetElt<Cplx,Index>& s, const ITReal& d, ManageStore& m);

void
doTask(const FillReal& f, ITReal& d);

void
doTask(const FillCplx& f, const ITReal& d, ManageStore& m);

void
doTask(const MultCplx& M, const ITReal& d, ManageStore& m);

void
doTask(const MultReal& m, ITReal& d);

Real
doTask(NormNoScale, const ITReal& d);

void
doTask(Conj,const ITReal& d);

void
doTask(TakeReal, const ITReal& );

void
doTask(TakeImag, const ITReal& d, ManageStore& m);

void
doTask(PrintIT<Index>& P, const ITReal& d);

Cplx
doTask(SumEls<Index>, const ITReal& d);

void
doTask(Write& W, const ITReal& d);

void
doTask(Contract<Index>& C,
       const ITReal& a1,
       const ITReal& a2,
       ManageStore& m);

void
doTask(const PlusEQ<Index>& P,
       ITReal& a1,
       const ITReal& a2);

bool inline
doTask(CheckComplex, const ITReal& d) { return false; }

} //namespace itensor

#endif

