//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CONDENSER_H
#define __ITENSOR_CONDENSER_H
#include "itensor/iqtensor.h"

namespace itensor {

struct IndexMap;

//
// Condenser
//
// Within one IQIndex, combine all Index's having the same QNs
//

class Condenser
    {
    public:

    Condenser() { }

    Condenser(const IQIndex& bigindex, IQIndex& smallindex);

    Condenser(const IQIndex& bigindex, const std::string& smallind_name);

    void 
    product(const IQTensor& t, IQTensor& res) const;

    IQTensor 
    operator*(const IQTensor& t) { IQTensor res; product(t,res); return res; }

    void 
    dag();

    void 
    prime(int inc = 1) { prime(All,inc); }

    void 
    prime(IndexType type = All, int inc = 1);

    const IQIndex& 
    bigind() const   { return bigind_; }

    const IQIndex& 
    smallind() const { return smallind_; }

    const std::vector<IndexMap>&
    maps() const { return maps_; }

    private:

    ///////////////
    IQIndex bigind_,   //uncondensed
            smallind_; //condensed

    std::vector<IndexMap> maps_;
    //////////////

    void 
    init(const std::string& smallind_name);

    }; //class Condenser

IQTensor inline
operator*(const IQTensor& t, const Condenser& c) 
    { 
    IQTensor res; 
    c.product(t,res); 
    return res; 
    }

Condenser inline
dag(Condenser res) { res.dag(); return res; }

std::ostream& 
operator<<(std::ostream & s, const Condenser & c);

struct IndexMap
    {
    IndexMap(const Index& small_,
             int i_,
             const Index& big_)
        :
        small(small_),
        i(i_),
        big(big_)
        { }

    Index small;
    int i;
    Index big;
    };

} //namespace itensor

#endif
