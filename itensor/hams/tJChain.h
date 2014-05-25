//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_TJCHAIN_H
#define __ITENSOR_HAMS_TJCHAIN_H
#include "../mpo.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

namespace itensor {

class tJChain 
    {
    public:

    tJChain(const SiteSet& sites,
            const OptSet& opts = Global::opts());

    Real
    t() const { return t_; }

    void
    t(Real val) { initted = false; t_ = val; }

    Real
    J() const { return J_; }

    void
    J(Real val) { initted = false; J_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    const SiteSet& sites_;
    Real t_,J_;
    bool initted;
    MPO H;

    void init_();

    }; //class tJChain

inline tJChain::
tJChain(const SiteSet& sites, const OptSet& opts)
    : 
    sites_(sites), 
    initted(false)
    { 
    J_ = opts.getReal("J",0.35);
    t_ = opts.getReal("t",1);
    }


void inline tJChain::
init_()
    {
    if(initted) return;
    Cout << "J is " << J_ << Endl;

    H = MPO(sites_);

    const int Ns = sites_.N();
    const int k = 10;

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("tjl",l),k);

    ITensor W;
    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Anc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(sites_.si(n),sites_.siP(n),row,col);

	// fermiPhase will be needed for longer range hopping, but not for this nn chain

        //W += sites_.Nupdn(n) * row(k) * col(1) * U_;	
        //W += multSiteOps(sites_.fermiPhase(n),sites_.Cup(n)) * row(k) * col(2) * t_;
        //W += multSiteOps(sites_.fermiPhase(n),sites_.Cdn(n)) * row(k) * col(3) * t_;
        //W += multSiteOps(sites_.Cdagup(n),sites_.fermiPhase(n)) * row(k) * col(4) * t_;
        //W += multSiteOps(sites_.Cdagdn(n),sites.fermiPhase(n)) * row(k) * col(5) * t_;

        W += sites_.Cup(n) * row(k) * col(2) * t_;
        W += sites_.Cdn(n) * row(k) * col(3) * t_;
        W += sites_.Cdagup(n) * row(k) * col(4) * t_;
        W += sites_.Cdagdn(n) * row(k) * col(5) * t_;

        W += sites_.sz(n) * row(k) * col(6) * J_;
        W += sites_.sp(n) * row(k) * col(7) * J_/2;
        W += sites_.sm(n) * row(k) * col(8) * J_/2;
        W += sites_.Ntot(n) * row(k) * col(9) * (-J_/4);
        W += sites_.id(n) * row(k) * col(k);

        W += sites_.id(n) * row(1) * col(1);
        W += sites_.Cdagup(n) * row(2) * col(1) * (-1.0);
        W += sites_.Cdagdn(n) * row(3) * col(1) * (-1.0);
        W += sites_.Cup(n) * row(4) * col(1) * (-1.0);
        W += sites_.Cdn(n) * row(5) * col(1) * (-1.0);
        W += sites_.sz(n) * row(6) * col(1);
        W += sites_.sm(n) * row(7) * col(1);
        W += sites_.sp(n) * row(8) * col(1);
        W += sites_.Ntot(n) * row(9) * col(1);
        }

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted = true;
    }

}; //namespace itensor

#undef Cout
#undef Endl
#undef Format

#endif
