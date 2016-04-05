//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QMIXED_H
#define __ITENSOR_QMIXED_H

#include "itensor/itdata/task_types.h"
#include "itensor/itdata/itdata.h"
#include "itensor/iqindex.h"

namespace itensor {

//template<typename index_type> 
//class ITensorT;
//using ITensor  = ITensorT<Index>;

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

const char*
typeNameOf(QMixed<Real> const& d);
const char*
typeNameOf(QMixed<Cplx> const& d);

template<typename V>
Cplx 
doTask(GetElt<IQIndex> const& g, QMixed<V> const& d);

template<typename E, typename T>
void
doTask(SetElt<E,IQIndex> const& S, QMixed<T> const& d, ManageStore & m);

//implementation of doTask(ToITensor...) is in iqtensor.cc
template<typename V>
ITensor
doTask(ToITensor & T, QMixed<V> const& d);

template<typename T>
void
doTask(PrintIT<IQIndex>& P, QMixed<T> const& D);


} //namespace itensor

#endif
