//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITREAL_H
#define __ITENSOR_ITREAL_H

#include "itensor/itdata/itdata.h"

namespace itensor {

//
// Optimization TODO: 
//  replace std::vector storage with
//  storage type only holding data ptr
//  and size
//

struct ITReal
    {
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


} //namespace itensor

#endif

