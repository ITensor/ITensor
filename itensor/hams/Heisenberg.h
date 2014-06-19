//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_HEISENBERG_H
#define __ITENSOR_HAMS_HEISENBERG_H

#include "../mpo.h"

namespace itensor {

class Heisenberg
    {
    public:

    Heisenberg(const SiteSet& sites, 
               const OptSet& opts = Global::opts());

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
         open_;
    IQMPO H;

    //
    //////////////////

    void 
    init_();

    }; //class Heisenberg

inline Heisenberg::
Heisenberg(const SiteSet& sites, const OptSet& opts)
    : 
    sites_(sites), 
    initted_(false)
    { 
    Ny_ = opts.getInt("Ny",1);
    Nx_ = sites_.N()/Ny_;
    J_ = opts.getReal("J",1.);
    Boundary_h_ = opts.getReal("Boundary_h",0.);
    Jz_ = opts.getReal("Jz",J_);
    infinite_ = opts.getBool("Infinite",false);
    open_ = opts.getBool("Open",false);
    }


void inline Heisenberg::
init_()
    {
    if(initted_) return;

    H = IQMPO(sites_);

    const int Ns = sites_.N(),
              kpm = Ny_,
              kd = Ny_+2,
              //k = kd+2*kpm,
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
                              qM[l],QN(+2),Out);
        }

    const IQIndex& last = (infinite_ ? links.at(0) : links.at(Ns));

    for(int n = 1; n <= Ns; ++n)
        {
        IQTensor& W = H.Anc(n);
        const
        IQIndex &row = links.at(n-1), 
                &col = (n==Ns ? last : links.at(n));

        W = IQTensor(dag(sites_.si(n)),sites_.siP(n),dag(row),col);

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
        const int x = (n-1)/Ny_+1, y = (n-1)%Ny_+1;
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
        if(!open_ && y == 1 && Ny_ >= 3)
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

    const IQTensor LH(links.at(0)(start)),
                   RH(dag(last)(1));

    if(infinite_)
        {
        H.Anc(0) = LH;
        H.Anc(Ns+1) = RH;
        }
    else
        {
        H.Anc(1) *= LH;
        H.Anc(Ns) *= RH;
        }

    initted_ = true;
    }

/*
//
//Older non-quantum-number version
//
void inline Heisenberg::
init_()
    {
    if(initted_) return;

    H = MPO(sites_);

    const int Ns = sites_.N();
    const int nop = 3;
    const int max_mpo_dist = Ny_;
    const int k = 5+nop*(max_mpo_dist-1);

    std::vector<Index> links(sites_.N()+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Anc(n);
        Index &row = links.at(n-1), &col = links.at(n);

        W = ITensor(sites_.si(n),sites_.siP(n),row,col);

        W += sites_.id(n) * row(1) * col(1);
        W += sites_.id(n) * row(k) * col(k);

        W += sites_.sz(n) * row(2) * col(1);
        W += sites_.sp(n) * row(3) * col(1);
        W += sites_.sm(n) * row(4) * col(1);

        //Horizontal bonds
        int mpo_dist = Ny_; 
        W += sites_.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * Jz_;
        W += sites_.sm(n) * row(k) * col(3+nop*(mpo_dist-1)) * J_/2;
        W += sites_.sp(n) * row(k) * col(4+nop*(mpo_dist-1)) * J_/2;

        //Add boundary field if requested
        const int x = (n-1)/Ny_+1, y = (n-1)%Ny_+1;
        if(Boundary_h_ != 0 && (x == 1 || x == Nx_))
            {
            Real eff_h = Boundary_h_;
            if(J_ > 0) eff_h *= (x%2==1 ? -1 : 1)*(y%2==1 ? -1 : 1);
            //cerr << format("Doing a staggered bf of %.2f at site %d (%d,%d)\n")%eff_h%n%x%y;
            W += sites_.sz(n) * ITensor(row(k),col(1)) * eff_h;
            }

        //The following only apply if Ny_ > 1:

        //String of identity ops
        for(int q = 1; q <= nop*(max_mpo_dist-1); ++q)
            { W += sites_.id(n) * row(1+nop+q) * col(1+q); }

        //Periodic BC bond (only for width 3 ladders or greater)
        if(y == 1 && Ny_ >= 3)
            {
            int mpo_dist = Ny_-1; 
            W += sites_.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * J_;
            W += sites_.sm(n) * row(k) * col(3+nop*(mpo_dist-1)) * J_/2;
            W += sites_.sp(n) * row(k) * col(4+nop*(mpo_dist-1)) * J_/2;
            }

        //N.N. bond along column
        if(y != Ny_)
            {
            W += sites_.sz(n) * row(k) * col(2) * J_;
            W += sites_.sm(n) * row(k) * col(3) * J_/2;
            W += sites_.sp(n) * row(k) * col(4) * J_/2;
            }
        }

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }
*/

}; //namespace itensor

#endif
