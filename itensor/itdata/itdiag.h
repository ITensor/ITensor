//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDIAG_H
#define __ITENSOR_ITDIAG_H

#include "itensor/itdata/itreal.h"
#include "itensor/detail/call_rewrite.h"

namespace itensor {

template<typename T>
class ITDiag : public RegisterData<ITDiag<T>>
    {
    public:
    using storage_type = typename std::vector<T>;
    using size_type = typename storage_type::size_type;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;
    using value_type = T;

    //
    // Data members
    //

    T val = 0;
    storage_type store;

    //
    // Constructors
    //

    ITDiag() { }

    template<typename InputIterator>
    ITDiag(InputIterator b, InputIterator e)
        :
        store(b,e)
        { }

    ITDiag(size_t size, T val)
        :
        store(size,val)
        { }

    ITDiag(T t) 
        : val(t) 
        { }

    explicit
    ITDiag(storage_type&& data)
        :
        store(std::move(data))
        { }

    virtual
    ~ITDiag() { }

    //
    // Accessors
    //

    bool
    allSame() const { return store.empty(); }

    //
    // std container like methods
    //

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

template<typename T>
void
read(std::istream& s, ITDiag<T>& dat)
    {
    read(s,dat.val);
    read(s,dat.store);
    }

template<typename T>
void
write(std::ostream& s, const ITDiag<T>& dat)
    {
    write(s,dat.val);
    write(s,dat.store);
    }

template <typename F, typename T,
          typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
void
doTask(ApplyIT<F>& A, ITDiag<T>& d) 
    { 
    if(d.allSame()) 
        {
        d.val = detail::call<T>(A.f,d.val);
        }
    else
        {
        for(auto& elt : d.store) elt = detail::call<T>(A.f,elt);
        }
    }

template <typename T>
Cplx
doTask(const GetElt<Index>& g, const ITDiag<T>& d);

void
doTask(Contract<Index>& C,
       const ITReal& t,
       const ITDiag<Real>& d,
       ManagePtr& mp);
void
doTask(Contract<Index>& C,
       const ITDiag<Real>& d,
       const ITReal& t,
       ManagePtr& mp);

void
doTask(const PlusEQ<Index>& P,
       ITDiag<Real>& a1,
       const ITDiag<Real>& a2);

template<typename T>
void
doTask(const FillReal& f, const ITDiag<T>& d, ManagePtr& mp);

template<typename T>
void
doTask(const FillCplx& f, const ITDiag<T>& d, ManagePtr& mp);

template<typename T>
void
doTask(const MultReal& m, ITDiag<T>& d);

template<typename T>
Real
doTask(const NormNoScale<Index>& N, const ITDiag<T>& d);

void
doTask(Conj, ITDiag<Cplx>& d);

void
doTask(Conj,const ITDiag<Real>& d);

void
doTask(TakeReal,const ITDiag<Cplx>& d, ManagePtr& mp);

void
doTask(TakeReal, const ITDiag<Real>& );

void
doTask(TakeImag,const ITDiag<Cplx>& d, ManagePtr& mp);

template<typename T>
void
doTask(PrintIT<Index>& P, const ITDiag<T>& d);

bool
doTask(CheckComplex,const ITDiag<Real>& d);

bool
doTask(CheckComplex,const ITDiag<Cplx>& d);

template <class T>
Cplx
doTask(SumEls<Index> S, const ITDiag<T>& d);

void
doTask(Write& W, const ITDiag<Real>& d);

void
doTask(Write& W, const ITDiag<Cplx>& d);

} //namespace itensor

#endif

