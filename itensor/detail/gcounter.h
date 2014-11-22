//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#pragma once

#include "autovector.h"

namespace itensor {
namespace detail {

class GCounter	// General Counter
    {
    public:

    // first are initial values of each index
    // last are final values of each index
    // i are the current values of each index
    autovector<long> first, 
                     last, 
                     i;

    // position of this GCounter value; starts at 0
    long ind;

    // for a GCounter that has indices i[1] through i[8], and starts its counts at 0, 
    // firstind = 1, lastind = 8, firstval = 0

    GCounter(long firstind, 
             long lastind, 
             long firstval = 1) 
        : 
        first(firstind,lastind,firstval), 
        last(firstind,lastind,firstval),
        i(firstind,lastind,firstval), 
        ind(0)
        { }

    // After constructing a GCounter g, calling g.setInd(j,s,e)
    // lets g.i[j] = s,s+1,...,e when iterating g
    void 
    setInd(long j, long f, long l)
        {
        first.ref(j) = f;
        last.ref(j) = l;
        i.ref(j) = f;
        ind = 0;
        }

    void 
    reset()
        {
        i = first;
        ind = 0;
        }

    GCounter& 
    operator++()
        {
        long mi = first.mini(), 
             ma = first.maxi();
        ++i.fastref(mi);
        ++ind;
        if(i.fast(mi) > last.fast(mi))
            {
            for(int j = mi+1; j <= ma; ++j)
                {
                i.fastref(j-1) = first.fast(j-1);
                ++i.fastref(j);
                if(i.fast(j) <= last.fast(j)) return *this;
                }
            i.fastref(mi) = first.fast(mi) - 1;	  // all done if get here; set !notdone()
            }
        return *this;
        }

    bool 
    notDone()
        {
        return i.fast(first.mini()) >= first.fast(first.mini());
        }
    };


};
};
