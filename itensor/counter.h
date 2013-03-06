//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_COUNTER_H
#define __ITENSOR_COUNTER_H

#include "indexset.h"

#define Array boost::array
#define Cout std::cout
#define Endl std::endl

//
// Counter
//

//
// A Counter "C" is iterated as
// for(; C.notDone(); ++C) { ... }
//
// When iterated, each index C.i[j]
// runs from 1 to C.n[j] such that taken
// together the C.i run over all values
// compatible with their ranges.
//
// The index C.i[1] changes every step,
// C.i[2] after every C.n[1] steps, etc.
//
// For example: (if C.n[1] == C.n[2] == 2, say)
// C.i[1]  C.i[2]   C.i[3]
//   1       1        1
//   2       1        1
//   1       2        1
//   2       2        1
//   1       1        2
//   2       1        2
//   ...
// 
// Counters also keep a straight count of which
// step they are on through the field C.ind
//
// ITensor data is ordered and the method
// ITensor::_ind is defined such that
// v(C.ind) == v(_ind(C.i[1],C.i[2],...,C.i[8])
// where v is the Vector in an ITDat.
//

class Counter
    {
    public:

    Array<int,NMAX+1> n, i;
    int ind;
    int rn_,r_;

    Counter();

    template <class IndexT>
    Counter(const Array<IndexT,NMAX+1>& ii,int rn,int r);

    template <class IndexT> 
    explicit
    Counter(const IndexSet<IndexT>& is);

    Counter& 
    operator++();

    //bool 
    //operator!=(const Counter& other) const;

    //bool 
    //operator==(const Counter& other) const;

    bool 
    notDone() const { return i[1] != 0; }

    void 
    reset(int a);

    };


//
// Counter implementation
//

inline 
Counter::
Counter() : rn_(0)
    {
    n.assign(1); n[0] = 0;
    reset(0);
    }

void inline Counter::
reset(int a)
    {
    i.assign(a);
    ind = 1;
    }

template <class IndexT>
Counter::
Counter(const Array<IndexT,NMAX+1>& ii, int rn, int r)
    {
    rn_ = rn;
    r_ = r;
    n[0] = 0;
    for(int j = 1; j <= rn_; ++j) 
        n[j] = ii[j].m();
    for(int j = rn_+1; j <= NMAX; ++j) 
        n[j] = 1;
    reset(1);
    }

template <class IndexT>
Counter::
Counter(const IndexSet<IndexT>& is)
    {
    rn_ = is.rn();
    r_ = is.r();
    n[0] = 0;
    for(int j = 1; j <= rn_; ++j) 
        n[j] = is.index(j).m();
    for(int j = rn_+1; j <= NMAX; ++j) 
        n[j] = 1;
    reset(1);
    }

inline
Counter& Counter::
operator++()
    {
    ++ind;
    ++i[1];
    if(i[1] > n[1])
    for(int j = 2; j <= rn_; ++j)
        {
        i[j-1] = 1;
        ++i[j];
        if(i[j] <= n[j]) break;
        }
    //set 'done' condition
    if(i[rn_] > n[rn_]) reset(0);
    return *this;
    }

/*
bool inline Counter::
operator!=(const Counter& other) const
    {
    for(int j = 1; j <= NMAX; ++j)
        { if(i[j] != other.i[j]) return true; }
    return false;
    }

bool inline Counter::
operator==(const Counter& other) const
    { return !(*this != other); }
    */

inline
std::ostream&
operator<<(std::ostream& s, const Counter& c)
    {
    s << "("; 
    for(int i = 1; i < c.r_; ++i)
        s << c.i[i] << " ";
    s << c.i[c.r_] << ")";
    return s;
    }

#undef Array
#undef Cout
#undef Endl
#endif
