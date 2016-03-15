//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QMIXED_H
#define __ITENSOR_QMIXED_H

#include "itensor/itdata/task_types.h"
#include "itensor/itdata/itdata.h"

namespace itensor {

template<typename T>
class QMixed
    {
    public:
    using value_type = T;
    using storage_type = std::vector<value_type>;
    using size_type = typename storage_type::size_type;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    storage_type store;

    QMixed() { }

    explicit
    QMixed(size_t size) : store(size) { }

    template<typename InputIterator>
    QMixed(InputIterator b, InputIterator e) : store(b,e) { }

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

template<typename V>
Cplx 
doTask(GetElt<IQIndex> const& g, QMixed<V> const& d)
    {
    return d[offset(g.is,g.inds)];
    }

namespace detail {
    template<typename E, typename T>
    struct SetEltHelper
        {
        void static
        set(SetElt<E,IQIndex> const& S, QMixed<T> const& D, ManageStore& m)
            {
            auto& Dnc = *m.modifyData(D);
            Dnc[offset(S.is,S.inds)] = S.elt;
            }
        };
    template<>
    struct SetEltHelper<Cplx,Real>
        {
        void static
        set(SetElt<Cplx,IQIndex> const& S, QMixed<Real> const& D, ManageStore & m)
            {
            auto& nd = *m.makeNewData<QMixed<Cplx>>(D.begin(),D.end());
            nd[offset(S.is,S.inds)] = S.elt;
            }
        };
} //namespace detail

template<typename E, typename T>
void
doTask(SetElt<E,IQIndex> const& S, QMixed<T> const& d, ManageStore & m)
    {
    detail::SetEltHelper<E,T>::set(S,d,m);
    }

//implementation of doTask(ToITensor...) is in iqtensor.cc
template<typename V>
ITensor
doTask(ToITensor & T, QMixed<V> const& d);

template<typename T>
void
doTask(PrintIT<IQIndex>& P, QMixed<T> const& D)
    {
    auto name = std::is_same<T,Real>::value ? "QMixed Real"
                                            : "QMixed Cplx";
    P.printInfo(D,name);
     
    auto r = rank(P.is);
    if(r == 0) 
        {
        P.s << "  ";
        P.s << formatVal(P.scalefac*D.store.front()) << "\n";
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(r);
    for(auto i : range(r))
        gc.setRange(i,0,P.is.extent(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = P.scalefac*D[offset(P.is,gc.i)];
        if(std::norm(val) >= Global::printScale())
            {
            P.s << "(";
            for(auto ii : range1(gc.i.mini(),gc.i.maxi()))
                {
                P.s << (1+gc[ii]);
                if(ii < gc.i.maxi()) P.s << ",";
                }
            P.s << ") ";

            P.s << formatVal(val) << "\n";
            }
        }
    }


} //namespace itensor

#endif
