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
               Args const& args = Global::args());

    operator MPO() { init_(); return H.toMPO(); }

    operator IQMPO() { init_(); return H; }

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
    IQMPO H;

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
    N_ = sites_.N();
    J_ = args.getReal("J",1.);
    Jz_ = args.getReal("Jz",J_);
    infinite_ = args.getBool("Infinite",false);
    }

void inline Heisenberg::
init_()
    {
    if(initted_) return;

    H = IQMPO(sites_);

    std::vector<IQIndex> links(N_+1);

    //The names of these indices refer to their Sz quantum numbers
    std::vector<Index> q0(N_+1),
                       qP(N_+1),
                       qM(N_+1);

    for(int l = 0; l <= N_; ++l) 
        {
        q0.at(l) = Index(nameint("q0_",l),3);
        qP.at(l) = Index(nameint("qP_",l),1);
        qM.at(l) = Index(nameint("qM_",l),1);

        links.at(l) = IQIndex(nameint("hl",l),
                              q0[l],QN( 0),
                              qP[l],QN(-2),
                              qM[l],QN(+2),
                              Out);
        }

    IQIndex const& last = (infinite_ ? links.at(0) : links.at(N_));

    for(int n = 1; n <= N_; ++n)
        {
        auto& W = H.Aref(n);
        auto row = dag(links.at(n-1));
        auto col = (n==N_ ? last : links.at(n));

        W = IQTensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * row(1) * col(1); //ending state
        W += sites_.op("Id",n) * row(2) * col(2); //starting state

        W += sites_.op("Sz",n) * row(3) * col(1);
        W += sites_.op("Sz",n) * row(2) * col(3) * Jz_;

        W += sites_.op("Sm",n) * row(4) * col(1);
        W += sites_.op("Sp",n) * row(2) * col(4) * J_/2;

        W += sites_.op("Sp",n) * row(5) * col(1);
        W += sites_.op("Sm",n) * row(2) * col(5) * J_/2;
        }

    auto LH = setElt(links.at(0)(2));
    auto RH = setElt(dag(last)(1));

    if(not infinite_)
        {
        //Multiply first and last
        //MPO tensor by edge vectors
        H.Aref(1) *= LH;
        H.Aref(N_) *= RH;
        }
    else
        {
        //Store edge vectors just before
        //and after first and last sites
        H.Aref(0) = LH;
        H.Aref(N_+1) = RH;
        }

    initted_ = true;
    }

} //namespace itensor

#endif
