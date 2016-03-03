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

    const SiteSet& sites_;
    int Ny_,
        Nx_;
    Real J_, 
         Boundary_h_,
         Jz_;
    bool initted_,
         infinite_,
         open_y_;
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
    Ny_ = args.getInt("Ny",1);
    Nx_ = sites_.N()/Ny_;
    J_ = args.getReal("J",1.);
    Boundary_h_ = args.getReal("Boundary_h",0.);
    Jz_ = args.getReal("Jz",J_);
    infinite_ = args.getBool("Infinite",false);
    open_y_ = args.getBool("Open",false);
    }

void inline Heisenberg::
init_()
    {
    if(initted_) return;

    H = IQMPO(sites_);

    int Ns = sites_.N(),
        kpm = Ny_,
        kd = Ny_+2,
        start = 2;

    std::vector<IQIndex> links(Ns+1);

    //The names of these indices refer to their Sz quantum numbers
    std::vector<Index> q0(Ns+1),
                       qP(Ns+1),
                       qM(Ns+1);

    for(int l = 0; l <= Ns; ++l) 
        {
        q0.at(l) = Index(nameint("q0_",l),kd);
        qP.at(l) = Index(nameint("qP_",l),kpm);
        qM.at(l) = Index(nameint("qM_",l),kpm);

        links.at(l) = IQIndex(nameint("hl",l),
                              q0[l],QN( 0),
                              qP[l],QN(-2),
                              qM[l],QN(+2),
                              Out);
        }

    IQIndex const& last = (infinite_ ? links.at(0) : links.at(Ns));

    for(int n = 1; n <= Ns; ++n)
        {
        auto& W = H.Anc(n);
        auto row = dag(links.at(n-1));
        auto col = (n==Ns ? last : links.at(n));

        W = IQTensor(dag(sites_.si(n)),sites_.siP(n),row,col);

        W += sites_.op("Id",n) * row(1) * col(1); //ending state
        W += sites_.op("Id",n) * row(2) * col(2); //starting state

        //
        //Sz Sz terms
        //
        W += sites_.op("Sz",n) * row(3) * col(1);
        //Horizontal bonds
        W += sites_.op("Sz",n) * row(start) * col(kd) * Jz_;

        //
        //S+ S- terms
        //
        W += sites_.op("Sm",n) * row(kd+1) * col(1);
        //Horizontal bonds
        W += sites_.op("Sp",n) * row(start) * col(kd+kpm) * J_/2;

        //
        //S- S+ terms
        //
        W += sites_.op("Sp",n) * row(kd+kpm+1) * col(1);
        //Horizontal bonds
        W += sites_.op("Sm",n) * row(start) * col(kd+2*kpm) * J_/2;

        //Add boundary field if requested
        auto x = (n-1)/Ny_+1; 
        auto y = (n-1)%Ny_+1;
        if(Boundary_h_ != 0 && (x == 1 || x == Nx_))
            {
            Real eff_h = Boundary_h_;
            if(J_ > 0) eff_h *= (x%2==1 ? -1 : 1)*(y%2==1 ? -1 : 1);
            printfln("Applying staggered h of %.2f at site %d (%d,%d)",eff_h,n,x,y);
            W += sites_.op("Sz",n) * row(start) * col(1) * eff_h;
            }

        //The following only apply if Ny_ > 1:

        //Strings of off-diagonal identity ops
        for(int q = 1; q <= Ny_-1; ++q)
            { 
            W += sites_.op("Id",n) * row(3+q) * col(2+q); 
            W += sites_.op("Id",n) * row(kd+1+q) * col(kd+q); 
            W += sites_.op("Id",n) * row(kd+kpm+1+q) * col(kd+kpm+q); 
            }

        //Periodic BC bond (only for width 3 ladders or greater)
        if(!open_y_ && y == 1 && Ny_ >= 3)
            {
            W += sites_.op("Sz",n) * row(start) * col(kd-1) * Jz_;
            W += sites_.op("Sp",n) * row(start) * col(kd+kpm-1) * J_/2;
            W += sites_.op("Sm",n) * row(start) * col(kd+2*kpm-1) * J_/2;
            }

        //N.N. bond along column
        if(y != Ny_)
            {
            W += sites_.op("Sz",n) * row(start) * col(3) * Jz_;
            W += sites_.op("Sp",n) * row(start) * col(kd+1) * J_/2;
            W += sites_.op("Sm",n) * row(start) * col(kd+kpm+1) * J_/2;
            }
        }

    auto LH = pick(links.at(0)(start));
    auto RH = pick(dag(last)(1));

    if(not infinite_)
        {
        //Multiply first and last
        //MPO tensor by edge vectors
        H.Anc(1) *= LH;
        H.Anc(Ns) *= RH;
        }
    else
        {
        //Store edge vectors just before
        //and after first and last sites
        H.Anc(0) = LH;
        H.Anc(Ns+1) = RH;
        }

    initted_ = true;
    }

} //namespace itensor

#endif
