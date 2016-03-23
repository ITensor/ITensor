//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#define __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#include "itensor/mps/mpo.h"

namespace itensor {

class ExtendedHubbard
    {
    public:

    ExtendedHubbard(const SiteSet& sites, 
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

    operator IQMPO() { init_(); return H.toIQMPO(); }

    private:

    //////////////////
    //
    // Data Members

    const SiteSet& sites_;
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
ExtendedHubbard(const SiteSet& sites,
                const OptSet& opts)
    :
    sites_(sites), 
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

    H = MPO(sites_);

    const int Ns = sites_.N();
    const int k = 7 + (t2_ == 0 ? 0 : 4);

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;
    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Anc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(sites_.si(n),sites_.siP(n),row,col);

        //Identity strings
        W += sites_.op("Id",n) * row(1) * col(1);
        W += sites_.op("Id",n) * row(k) * col(k);

        //Hubbard U term
        W += sites_.op("Nupdn",n) * row(k) * col(1) * U_;

        //Hubbard V1 term
        W += sites_.op("Ntot",n) * row(k-1) * col(1);
        W += sites_.op("Ntot",n) * row(k) * col(k-1) * V1_;

        //Kinetic energy/hopping terms, defined as -t_*(c^d_i c_{i+1} + h.c.)
        W += sites_.op("Aup*F",n)    * row(k) * col(2) * t1_;
        W += sites_.op("Adn",n)      * row(k) * col(3) * t1_;
        W += sites_.op("Adagup*F",n) * row(k) * col(4) * -t1_;
        W += sites_.op("Adagdn",n)   * row(k) * col(5) * -t1_;

        if(t2_ != 0)
            {
            W += sites_.op("Aup*F",n)    * row(k) * col(6) * t2_;
            W += sites_.op("Adn",n)      * row(k) * col(7) * t2_;
            W += sites_.op("Adagup*F",n) * row(k) * col(8) * -t2_;
            W += sites_.op("Adagdn",n)   * row(k) * col(9) * -t2_;

            W += sites_.op("F",n)   * row(6) * col(2);
            W += sites_.op("F",n)   * row(7) * col(3);
            W += sites_.op("F",n)   * row(8) * col(4);
            W += sites_.op("F",n)   * row(9) * col(5);
            }

        W += sites_.op("Adagup",n)   * row(2) * col(1);
        W += sites_.op("F*Adagdn",n) * row(3) * col(1);
        W += sites_.op("Aup",n)      * row(4) * col(1);
        W += sites_.op("F*Adn",n)    * row(5) * col(1);

        }

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

} //namespace itensor

#endif
