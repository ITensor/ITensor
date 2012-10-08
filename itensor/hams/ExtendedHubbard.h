//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#define __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#include "../mpo.h"
#include "../hams.h"
#include "../model/hubbard.h"

class ExtendedHubbard : public MPOBuilder
    {
    public:
    typedef MPOBuilder Parent;

    ExtendedHubbard(const Hubbard& model, 
                    Real t = 1, Real U = 0, Real V = 0);

    Real
    t() const { return t_; }
    void
    t(Real val) { initted = false; t_ = val; }

    Real
    U() const { return U_; }
    void
    U(Real val) { initted = false; U_ = val; }

    Real
    V() const { return V_; }
    void
    V(Real val) { initted = false; V_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    const Hubbard& model_;
    Real t_,U_,V_;
    bool initted;
    MPO H;

    void init_();

    }; //class HubbardChain

inline HubbardChain::
HubbardChain(const Hubbard::Model& model, 
             Real t, Real U, Real V) 
    : Parent(model), 
      model_(model), 
      t_(t), 
      U_(U), 
      V_(V), 
      initted(false)
    { }

void inline HubbardChain::
init_()
    {
    if(initted) return;

    H = MPO(model_);

    const int k = 6;

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;
    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.AAnc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        W += model_.Nupdn(n) * row(k) * col(1) * U_;
        W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * t_;
        W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(3) * t_;
        W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(4) * t_;
        W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(5) * t_;
        W += model_.id(n) * row(k) * col(k);

        W += model_.id(n) * row(1) * col(1);
        W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
        W += model_.Cdagdn(n) * row(3) * col(1) * (-1.0);
        W += model_.Cup(n) * row(4) * col(1) * (-1.0);
        W += model_.Cdn(n) * row(5) * col(1) * (-1.0);
        }

    H.AAnc(1) *= makeLedge(links.at(0));
    H.AAnc(Ns) *= makeRedge(links.at(Ns));

    initted = true;
    }

#endif
