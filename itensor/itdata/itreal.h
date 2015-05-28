//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITREAL_H
#define __ITENSOR_ITREAL_H

#include "itensor/itdata/task_types.h"
#include "itensor/util/readwrite.h"

namespace itensor {

class ManagePtr;

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
    read(s,dat.store);
    }

void inline
write(std::ostream& s, const ITReal& dat)
    {
    write(s,dat.store);
    }


Cplx 
doTask(const GetElt<Index>& g, const ITReal& d);

void
doTask(const SetElt<Real,Index>& s, ITReal& d);

void
doTask(const SetElt<Cplx,Index>& s, const ITReal& d, ManagePtr& mp);

void
doTask(const FillReal& f, ITReal& d);

void
doTask(const FillCplx& f, const ITReal& d, ManagePtr& mp);

void
doTask(const MultCplx& M, const ITReal& d, ManagePtr& mp);

void
doTask(const MultReal& m, ITReal& d);

Real
doTask(const NormNoScale<Index>& N, const ITReal& d);

void
doTask(Conj,const ITReal& d);

void
doTask(TakeReal, const ITReal& );

void
doTask(TakeImag, const ITReal& d, ManagePtr& mp);

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
       ManagePtr& mp);

void
doTask(const PlusEQ<Index>& P,
       ITReal& a1,
       const ITReal& a2);

bool inline
doTask(CheckComplex, const ITReal& d) { return false; }

} //namespace itensor

#endif

