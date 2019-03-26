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
    MPO H_;
    public:

    ExpHeisenberg(SiteSet const& sites, 
             Cplx tau,
             Args const& args = Args::global());

    operator MPO() { init_(); return H_; }

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
    N_ = sites_.length();
    }

void inline ExpHeisenberg::
init_()
    {
    if(initted_) return;

    H_ = MPO(sites_);

    auto links = std::vector<Index>(N_+1);
    for(int l = 0; l <= N_; ++l) 
        {
        auto ts = format("Link,l=%d",l);
        links.at(l) = Index(QN({"Sz",0}),2,
                            QN({"Sz",-2}),1,
                            QN({"Sz",+2}),1,
                            ts);
        }

    for(int n = 1; n <= N_; ++n)
        {
        auto& W = H_.ref(n);
        auto row = dag(links.at(n-1));
        auto col = links.at(n);

        W = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * setElt(row(1)) * setElt(col(1));

        W += -tau_ * sites_.op("Sz",n) * setElt(row(1)) * setElt(col(2));
        W += sites_.op("Sz",n) * setElt(row(2)) * setElt(col(1));

        W += -(tau_/2.) * sites_.op("S+",n) * setElt(row(1)) * setElt(col(3));
        W += sites_.op("S-",n) * setElt(row(3)) * setElt(col(1));

        W += -(tau_/2.) * sites_.op("S-",n) * setElt(row(1)) * setElt(col(4));
        W += sites_.op("S+",n) * setElt(row(4)) * setElt(col(1));
        }

    H_.ref(1)  *= setElt(links.at(0)(1));
    H_.ref(N_) *= setElt(dag(links.at(N_)(1)));

    initted_ = true;
    }

} //namespace itensor

#endif
