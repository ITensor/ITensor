//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_HUBBARDCHAIN_H
#define __ITENSOR_HAMS_HUBBARDCHAIN_H
#include "../mpo.h"
#include "../hams.h"
#include "../model/hubbard.h"

class HubbardChain
    {
    public:

    HubbardChain(const Hubbard& model,
                 const OptSet& opts = Global::opts());

    Real
    t() const { return t_; }
    void
    t(Real val) { initted_ = false; t_ = val; }

    Real
    U() const { return U_; }
    void
    U(Real val) { initted_ = false; U_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    ///////////////////
    //
    // Data Members

    const Hubbard& model_;
    Real t_,U_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class HubbardChain

inline HubbardChain::
HubbardChain(const Hubbard& model, 
             const OptSet& opts)
    : 
    model_(model), 
    initted_(false)
    { 
    U_ = opts.getReal("U",0);
    t_ = opts.getReal("t",1);
    }

void inline HubbardChain::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.N();
    const int k = 6;

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;
    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Aref(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        //Identity strings
        W += model_.id(n) * row(1) * col(1);
        W += model_.id(n) * row(k) * col(k);

        //Hubbard U
        W += model_.Nupdn(n) * row(k) * col(1) * U_;

        //Kinetic energy/hopping terms, defined as -t_*(c^d_i c_{i+1} + h.c.)
        W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * t_;
        W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(3) * t_;
        W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(4) * t_;
        W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(5) * t_;
        W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
        W += model_.Cdagdn(n) * row(3) * col(1) * (-1.0);
        W += model_.Cup(n) * row(4) * col(1) * (-1.0);
        W += model_.Cdn(n) * row(5) * col(1) * (-1.0);
        }

    H.Aref(1) *= ITensor(links.at(0)(k));
    H.Aref(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

#endif
