//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LATTICEBOND_H__
#define __ITENSOR_LATTICEBOND_H__

#include <vector>
#include "itensor/global.h"

namespace itensor {

struct LatticeBond;

using LatticeGraph = std::vector<LatticeBond>;

struct LatticeBond
    {
    int s1 = 0,
        s2 = 0;
    std::string type;
    Real x1 = NAN,
         y1 = NAN,
         x2 = NAN,
         y2 = NAN;

    LatticeBond() { }

    LatticeBond(int s1_, int s2_)
      : s1{s1_}, 
        s2{s2_} 
        { }

    LatticeBond(int s1_, int s2_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_)
      : s1{s1_}, 
        s2{s2_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_}
        { }

    LatticeBond(int s1_, int s2_, std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        type{type_} 
        { }

    LatticeBond(int s1_, int s2_, 
                Real x1_,  Real y1_,
                Real x2_,  Real y2_,
                std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        type{type_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_}
        { }
    };

inline std::ostream& 
operator<<(std::ostream & s, LatticeBond const& b) 
    { 
    //s << format("(%*d,%*d",3,b.s1,3,b.s2);
    s << format("(%d,%d",b.s1,b.s2);
    if(b.type.size()!=0) s << "," << b.type;
    s << ")";
    if(!std::isnan(b.x1) && !std::isnan(b.y1))
        {
        s << format("[%s,%s",b.x1,b.y1);
        if(!std::isnan(b.x2) && !std::isnan(b.y2))
            {
            s << format(";%s,%s]",b.x2,b.y2);
            }
        else
            {
            s << "]";
            }
        }
    return s;
    }

inline std::ostream& 
operator<<(std::ostream& s, LatticeGraph const& G) 
    { 
    for(auto& b : G)
        {
        s << b << "\n";
        }
    return s;
    }

} //namespace itensor

#endif
