//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __EXP_HEISENBERG_H
#define __EXP_HEISENBERG_H

#include "itensor/mps/mpo.h"

namespace itensor {

class ExpHeisenberg
    {
    SiteSet const& sites_;
    int N_ = 0;
    Cplx tau_ = 0;
    bool initted_ = false;
    IQMPO H_;
    public:

    ExpHeisenberg(SiteSet const& sites, 
             Cplx tau,
             Args const& args = Global::args());

    operator IQMPO() { init_(); return H_; }

    private:

    void 
    init_();

    }; //class ExpHeisenberg

inline ExpHeisenberg::
ExpHeisenberg(SiteSet const& sites, 
         Cplx tau,
         Args const& args)
  : sites_(sites), 
    tau_(tau),
    initted_(false)
    { 
    N_ = sites_.N();
    }

void inline ExpHeisenberg::
init_()
    {
    if(initted_) return;

    H_ = IQMPO(sites_);

    auto links = std::vector<IQIndex>(N_+1);
    for(int l = 0; l <= N_; ++l) 
        {
        links.at(l) = IQIndex(nameint("hl",l),
                              Index(nameint("z",l),2),QN("Sz=",0),
                              Index(nameint("p",l),1),QN("Sz=",-2),
                              Index(nameint("m",l),1),QN("Sz=",+2));
        }

    for(int n = 1; n <= N_; ++n)
        {
        auto& W = H_.Aref(n);
        auto row = dag(links.at(n-1));
        auto col = links.at(n);

        W = IQTensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * row(1) * col(1);

        W += -tau_ * sites_.op("Sz",n) * row(1) * col(2);
        W += sites_.op("Sz",n) * row(2) * col(1);

        W += -(tau_/2.) * sites_.op("S+",n) * row(1) * col(3);
        W += sites_.op("S-",n) * row(3) * col(1);

        W += -(tau_/2.) * sites_.op("S-",n) * row(1) * col(4);
        W += sites_.op("S+",n) * row(4) * col(1);
        }

    H_.Aref(1)  *= setElt(links.at(0)(1));
    H_.Aref(N_) *= setElt(dag(links.at(N_)(1)));

    initted_ = true;
    }

} //namespace itensor

#endif
