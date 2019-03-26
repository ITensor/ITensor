//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_HEISENBERG_H
#define __ITENSOR_HAMS_HEISENBERG_H

#include "itensor/mps/mpo.h"

namespace itensor {

class Heisenberg
    {
    public:

    Heisenberg(SiteSet const& sites, 
               Args const& args = Args::global());

    operator MPO() { init_(); return H; }

    private:

    //////////////////
    //
    // Data Members

    SiteSet const& sites_;
    int N_;
    Real J_, 
         Jz_;
    bool initted_,
         infinite_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class Heisenberg

inline Heisenberg::
Heisenberg(SiteSet const& sites, 
           Args const& args)
  : sites_(sites), 
    initted_(false)
    { 
    N_ = sites_.length();
    J_ = args.getReal("J",1.);
    Jz_ = args.getReal("Jz",J_);
    infinite_ = args.getBool("Infinite",false);
    }

void inline Heisenberg::
init_()
    {
    if(initted_) return;

    H = MPO(sites_);

    std::vector<Index> links(N_+1);

    for(int l = 0; l <= N_; ++l) 
        {
        auto ts = format("Link,l=%d",l);
        links.at(l) = Index(QN({"Sz", 0}),3,
                            QN({"Sz",-2}),1,
                            QN({"Sz",+2}),1,
                            Out,
                            ts);
        }

    Index const& last = (infinite_ ? links.at(0) : links.at(N_));

    for(int n = 1; n <= N_; ++n)
        {
        auto& W = H.ref(n);
        auto row = dag(links.at(n-1));
        auto col = (n==N_ ? last : links.at(n));

        W = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * setElt(row(1)) * setElt(col(1)); //ending state
        W += sites_.op("Id",n) * setElt(row(2)) * setElt(col(2)); //starting state

        W += sites_.op("Sz",n) * setElt(row(3)) * setElt(col(1));
        W += sites_.op("Sz",n) * setElt(row(2)) * setElt(col(3)) * Jz_;

        W += sites_.op("Sm",n) * setElt(row(4)) * setElt(col(1));
        W += sites_.op("Sp",n) * setElt(row(2)) * setElt(col(4)) * J_/2;

        W += sites_.op("Sp",n) * setElt(row(5)) * setElt(col(1));
        W += sites_.op("Sm",n) * setElt(row(2)) * setElt(col(5)) * J_/2;
        }

    auto LH = setElt(links.at(0)(2));
    auto RH = setElt(dag(last)(1));

    if(not infinite_)
        {
        //Multiply first and last
        //MPO tensor by edge vectors
        H.ref(1) *= LH;
        H.ref(N_) *= RH;
        }
    else
        {
        //Store edge vectors just before
        //and after first and last sites
        H.ref(0) = LH;
        H.ref(N_+1) = RH;
        }

    initted_ = true;
    }

} //namespace itensor

#endif
