//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DENSE_H
#define __ITENSOR_DENSE_H

#include "itensor/itdata/task_types.h"
#include "itensor/util/readwrite.h"
#include "itensor/detail/call_rewrite.h"
#include "itensor/itdata/itdata.h"

namespace itensor {

template<typename T>
class Dense;

using DenseReal = Dense<Real>;
using DenseCplx = Dense<Cplx>;

template<typename T>
class Dense
    {
    static_assert(not std::is_const<T>::value,
                  "Template argument to Dense storage should not be const");
    public:
    using value_type = T;
    using storage_type = std::vector<value_type>;
    using size_type = typename storage_type::size_type;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    //
    // Data members
    //

    storage_type store;

    //
    // Constructors
    //

    Dense() { }

    explicit
    Dense(size_t size) : store(size) { }

    Dense(size_t size, value_type val) 
      : store(size,val)
        { }

    template<typename InputIterator>
    Dense(InputIterator b, InputIterator e) : store(b,e) { }

    Dense(storage_type&& data) : store(std::move(data)) { }

    //
    //std container like methods
    //

    value_type&
    operator[](size_type i) 
        { 
#ifdef DEBUG
        return store.at(i);
#else
        return store[i]; 
#endif
        }
    value_type const&
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

    value_type*
    data() { return store.data(); }
    value_type const*
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

const char*
typeNameOf(DenseReal const& d);
const char*
typeNameOf(DenseCplx const& d);

template<typename T>
bool constexpr
isReal(Dense<T> const& t) { return std::is_same<T,Real>::value; }

template<typename T>
bool constexpr
isCplx(Dense<T> const& t) { return std::is_same<T,Cplx>::value; }

Data inline
realData(DenseReal & d) { return Data(d.data(),d.size()); }

Datac inline
realData(DenseReal const& d) { return Datac(d.data(),d.size()); }

Data inline
realData(DenseCplx & d) { return Data(reinterpret_cast<Real*>(d.data()),2*d.size()); }

Datac inline
realData(DenseCplx const& d) { return Datac(reinterpret_cast<const Real*>(d.data()),2*d.size()); }

template<typename T>
void 
read(std::istream& s, Dense<T> & dat)
    {
    itensor::read(s,dat.store);
    }

template<typename T>
void
write(std::ostream& s, Dense<T> const& dat)
    {
    itensor::write(s,dat.store);
    }

template<typename F, typename T>
void
doTask(ApplyIT<F>& A, Dense<T> const& d, ManageStore & m)
    { 
    using new_type = ApplyIT_result_of<T,F>;
    if(switchesType<T>(A))
        {
        auto *nd = m.makeNewData<Dense<new_type>>(d.size());
        for(auto i : range(d))
            {
            A(d.store[i],nd->store[i]);
            }
        }
    else
        {
        auto *md = m.modifyData(d);
        for(auto& el : *md) A(el);
        }
    }

template<typename F, typename T>
void
doTask(VisitIT<F>& V, Dense<T> const& d)
    { 
    for(auto& elt : d) detail::call<void>(V.f,V.scale_fac * elt);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, DenseReal & D)
    {
    stdx::generate(D,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, DenseCplx const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<DenseReal>(D.size());
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, Dense<Real> const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<DenseCplx>(D.size());
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, Dense<Cplx> & D)
    {
    stdx::generate(D,G.f);
    }


Cplx 
doTask(GetElt<Index> const& g, DenseReal const& d);
Cplx 
doTask(GetElt<Index> const& g, DenseCplx const& d);

template<typename E, typename T>
void
doTask(SetElt<E,Index> const& S, Dense<T> const& d, ManageStore & m);

void
doTask(Fill<Real> const& f, DenseReal & d);

void
doTask(Fill<Real> const& f, DenseCplx const& d, ManageStore & m);

void
doTask(Fill<Cplx> const& f, DenseReal const& D, ManageStore & m);

void
doTask(Fill<Cplx> const& f, DenseCplx & D);

template<typename T>
void
doTask(Mult<Real> const& f, Dense<T> & D);

void
doTask(Mult<Cplx> const& f, Dense<Real> const& D, ManageStore & m);

void
doTask(Mult<Cplx> const& f, Dense<Cplx> & D);

template<typename T>
Real
doTask(NormNoScale, Dense<T> const& d);

void
doTask(Conj,DenseReal const& d);

void
doTask(Conj,DenseCplx & d);

void
doTask(TakeReal, DenseReal const& );

void
doTask(TakeReal, DenseCplx const& D, ManageStore & m);

void
doTask(TakeImag, DenseReal & d);

void
doTask(TakeImag, DenseCplx const& d, ManageStore & m);

void
doTask(MakeCplx, DenseReal const& d, ManageStore & m);

void inline
doTask(MakeCplx, DenseCplx const& d) { }

template<typename T>
bool constexpr
doTask(CheckComplex, Dense<T> const& d) { return isCplx(d); }

template<typename T>
void
doTask(PrintIT<Index>& P, Dense<T> const& d);

template<typename T>
Cplx
doTask(SumEls<Index>, Dense<T> const& d);

auto constexpr inline
doTask(StorageType const& S, DenseReal const& d) ->StorageType::Type { return StorageType::DenseReal; }

auto constexpr inline
doTask(StorageType const& S, DenseCplx const& d) ->StorageType::Type { return StorageType::DenseCplx; }

template<typename T1,typename T2>
void
doTask(Contract<Index> & C,
       Dense<T1> const& L,
       Dense<T2> const& R,
       ManageStore & m);

template<typename T1, typename T2>
void
doTask(NCProd<Index>& NCP,
       Dense<T1> const& D1,
       Dense<T2> const& D2,
       ManageStore& m);

template<typename T1, typename T2>
void
doTask(PlusEQ<Index> const& P,
       Dense<T1> const& D1,
       Dense<T2> const& D2,
       ManageStore & m);

template<typename T>
void
doTask(Order<Index> const& P,
       Dense<T> & dA);

template<typename T>
void
permuteDense(Permutation const& P,
             Dense<T>   const& dA,
             IndexSet   const& Ais,
             Dense<T>        & dB,
             IndexSet   const& Bis);

} //namespace itensor

#endif
