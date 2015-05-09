//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_COUNTER_H
#define __ITENSOR_COUNTER_H

#include "indexset.h"

#define Cout std::cout
#define Endl std::endl

namespace itensor {

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
//   0       0        0
//   1       0        0
//   0       1        0
//   1       1        0
//   0       0        1
//   1       0        1
//   ...
// 
// Counters also keep a straight count of which
// step they are on through the field C.ind
// (which starts from 0).
//
// ITensor data is ordered and the method
// ITensor::_ind is defined such that
// v[C.ind] == v[_ind(C.i[1],C.i[2],...,C.i[8]]
// where v is the Vector in an ITDat.
//

class Counter
    {
    public:

    array<int,NMAX+1> n, 
                      i;
    int ind,
        rn,
        r;

    Counter();

    template <class IndexT>
    Counter(const array<IndexT,NMAX+1>& ii,int rn,int r);

    template <class IndexT> 
    explicit
    Counter(const IndexSet<IndexT>& is);

    Counter& 
    operator++();

    bool 
    notDone() const { return i[1] >= 0; }

    void 
    reset();

    };


//
// Counter implementation
//

inline Counter::
Counter() 
    : rn(0), r(0)
    {
    n.fill(1); 
    n[0] = 0;
    reset();
    }

void inline Counter::
reset()
    {
    i.fill(0);
    ind = 0;
    }

template <class IndexT>
Counter::
Counter(const array<IndexT,NMAX+1>& ii, int rn_, int r_)
    {
    rn = rn_;
    r = r_;
    n[0] = 0;
    for(int j = 1; j <= rn; ++j) 
        n[j] = ii[j].m();
    for(int j = rn+1; j <= NMAX; ++j) 
        n[j] = 1;
    reset();
    }

template <class IndexT>
Counter::
Counter(const IndexSet<IndexT>& is)
    {
    rn = is.rn();
    r = is.r();
    n[0] = 0;
    for(int j = 1; j <= rn; ++j) 
        n[j] = is.index(j).m();
    for(int j = rn+1; j <= NMAX; ++j) 
        n[j] = 1;
    reset();
    }

inline
Counter& Counter::
operator++()
    {
    ++ind;
    ++i[1];
    if(i[1] >= n[1])
        {
        for(int j = 2; j <= rn; ++j)
            {
            i[j-1] = 0;
            ++i[j];
            if(i[j] < n[j]) break;
            }
        }
    //set 'done' condition
    if(i[rn] >= n[rn]) 
        {
        i[1] = -1;
        }
    return *this;
    }

inline
std::ostream&
operator<<(std::ostream& s, const Counter& c)
    {
    s << "("; 
    if(c.r >= 1) s << (c.i[1]+1);
    for(int i = 2; i <= c.r; ++i)
        {
        s << "," << (c.i[i]+1);
        }
    s << ")";
    return s;
    }

} //namespace itensor

#undef Cout
#undef Endl
#endif
