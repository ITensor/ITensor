//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITCPLX_H
#define __ITENSOR_ITCPLX_H

#include "itensor/itdata/itreal.h"
#include "itensor/util/safe_ptr.h"

namespace itensor {

//
// Optimization TODO: 
//  replace std::vector storage with
//  storage type only holding data ptr
//  and size
//

struct ITCplx
    {
    using storage_type = std::vector<Real>;
    using size_type = storage_type::size_type;
    using iterator = storage_type::iterator;
    using const_iterator = storage_type::const_iterator;

    //
    // Data members
    //

    storage_type store;

    //
    // Constructors
    //

    ITCplx() { }

    ITCplx(size_t size) : store(2*size,0) { }

    //Set all elements equal to (val,0)
    ITCplx(size_t size, Real val) 
        : 
        store(2*size,0) 
        { 
        std::fill(store.begin(),store.begin()+csize(),val);
        }
    //Set all elements equal to (val.real(),val.imag())
    ITCplx(size_t size, const Complex& val) 
        : 
        store(2*size) 
        { 
        fill(val);
        }

    explicit
    ITCplx(const ITReal& d)
        : 
        store(2*d.size(),0) 
        { 
        std::copy(d.begin(),d.end(),store.begin());
        }

    template<typename InputIterator>
    ITCplx(InputIterator b, InputIterator e) : store(b,e) { }

    //
    // Accessors
    //

    Complex
    get(size_type i) const
        {
        return Complex(store[i],store[csize()+i]);
        }

    void
    set(size_type i, const Complex& z) 
        {
        store[i] = z.real();
        store[csize()+i] = z.imag();
        }

    ITCplx&
    operator*=(const Complex& z)
        {
        auto r = MAKE_SAFE_PTR(rstart(),csize());
        auto re = istart(); //imag start is real end
        auto i = MAKE_SAFE_PTR(istart(),csize());
        auto a = z.real(),
             b = z.imag();
        for(; r != re; ++r, ++i)
            {
            auto nr = *r*a-*i*b;
            auto ni = *i*a+*r*b;
            *r = nr;
            *i = ni;
            }
        return *this;
        }

    void
    fill(const Complex& z)
        {
        std::fill(store.begin(),store.begin()+csize(),z.real());
        std::fill(store.begin()+csize(),store.end(),z.imag());
        }

    size_type
    csize() const { return store.size()/2; }

    Real*
    rstart() { return store.data(); }
    const Real*
    rstart() const { return store.data(); }
    const Real*
    rend() const { return store.data()+csize(); }

    Real*
    istart() { return store.data()+csize(); }
    const Real*
    istart() const { return store.data()+csize(); }
    const Real*
    iend() const { return store.data()+store.size(); }

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

void inline
read(std::istream& s, ITCplx& dat)
    {
    read(s,dat.store);
    }

void inline
write(std::ostream& s, const ITCplx& dat)
    {
    write(s,dat.store);
    }

} //namespace itensor

#endif

