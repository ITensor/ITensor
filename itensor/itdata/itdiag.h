//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDIAG_H
#define __ITENSOR_ITDIAG_H

#include "itdata.h"
#include "../simpletensor.h"

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

} //namespace itensor

#endif

