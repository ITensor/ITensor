//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#define __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#include "../mpo.h"
#include "../hams.h"
#include "../model/hubbard.h"

class ExtendedHubbard
    {
    public:

    ExtendedHubbard(const Hubbard& model, 
                    const OptSet& opts = Global::opts());

    Real
    t1() const { return t1_; }
    void
    t1(Real val) { initted_ = false; t1_ = val; }

    Real
    t2() const { return t2_; }
    void
    t2(Real val) { initted_ = false; t2_ = val; }

    Real
    U() const { return U_; }
    void
    U(Real val) { initted_ = false; U_ = val; }

    Real
    V1() const { return V1_; }
    void
    V1(Real val) { initted_ = false; V1_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    //////////////////
    //
    // Data Members

    const Hubbard& model_;
    Real U_,
         t1_,
         t2_,
         V1_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void init_();

    }; //class HubbardChain

inline ExtendedHubbard::
ExtendedHubbard(const Hubbard& model,
                const OptSet& opts)
    :
    model_(model), 
    initted_(false)
    { 
    U_ = opts.getReal("U",0);
    t1_ = opts.getReal("t1",1);
    t2_ = opts.getReal("t2",0);
    V1_ = opts.getReal("V1",0);
    }

void inline ExtendedHubbard::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.NN();
    const int k = 7 + (t2_ == 0 ? 0 : 4);

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;
    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.AAnc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        //Identity strings
        W += model_.id(n) * row(1) * col(1);
        W += model_.id(n) * row(k) * col(k);

        //Hubbard U term
        W += model_.Nupdn(n) * row(k) * col(1) * U_;

        //Hubbard V1 term
        W += model_.Ntot(n) * row(k-1) * col(1);
        W += model_.Ntot(n) * row(k) * col(k-1) * V1_;

        if(t2_ == 0)
            {
            //Kinetic energy/hopping terms, defined as -t_*(c^d_i c_{i+1} + h.c.)
            W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * t1_;
            W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(3) * t1_;
            W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(4) * t1_;
            W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(5) * t1_;

            W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
            W += model_.Cdagdn(n) * row(3) * col(1) * (-1.0);
            W += model_.Cup(n) * row(4) * col(1) * (-1.0);
            W += model_.Cdn(n) * row(5) * col(1) * (-1.0);
            }
        else // t2_ != 0
            {
            W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * t1_;
            W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(3) * t2_;
            W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(4) * t1_;
            W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(5) * t2_;
            W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(6) * t1_;
            W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(7) * t2_;
            W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(8) * t1_;
            W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(9) * t2_;

            W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);
            W += model_.fermiPhase(n)*row(3)*col(2);
            W += model_.Cdagdn(n)*row(4)*col(1)*(-1.0);
            W += model_.fermiPhase(n)*row(5)*col(4);
            W += model_.Cup(n)*row(6)*col(1)*(-1.0);
            W += model_.fermiPhase(n)*row(7)*col(6);
            W += model_.Cdn(n)*row(8)*col(1)*(-1.0);
            W += model_.fermiPhase(n)*row(9)*col(8);
            }
        }

    H.AAnc(1) *= ITensor(links.at(0)(k));
    H.AAnc(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

#endif
