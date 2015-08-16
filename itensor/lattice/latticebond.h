//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LATTICEBOND_H__
#define __ITENSOR_LATTICEBOND_H__

#include <vector>
#include "global.h"

namespace itensor {

struct LatticeBond;

using Lattice = std::vector<LatticeBond>;

struct LatticeBond
    {
    int s1 = 0,
        s2 = 0;
    std::string type;

    LatticeBond() { }

    LatticeBond(int s1_, int s2_)
      : s1{s1_}, 
        s2{s2_} 
        { }

    LatticeBond(int s1_, int s2_, std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        type{type_} 
        { }
    };

inline std::ostream& 
operator<<(std::ostream& s, const LatticeBond& b) 
    { 
    s << "(" << b.s1 << "," << b.s2;
    if(b.type.size()!=0) s << "," << b.type;
    s << ")";
    return s;
    }

} //namespace itensor

#endif
