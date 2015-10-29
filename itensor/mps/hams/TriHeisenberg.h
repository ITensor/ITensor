//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_TRIANGULARHEISENBERG_H
#define __ITENSOR_HAMS_TRIANGULARHEISENBERG_H
#include "itensor/mps/mpo.h"

namespace itensor {

class TriangularHeisenberg
    {
    public:

    TriangularHeisenberg(const SiteSet& sites, 
                         const OptSet& opts = Global::opts());

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H.toIQMPO(); }

    //------------------------------------------------------//

    private:

    const SiteSet& sites_;
    int Ny_,
        Nx_;
    Real Jz_, 
         Jxy_,
         boundary_h_;
    bool initted_;
    MPO H;

    void init_();

    }; //class TriangularHeisenberg

inline TriangularHeisenberg::
TriangularHeisenberg(const SiteSet& sites, const OptSet& opts)
    :
    sites_(sites), 
    initted_(false)
    { 
    Ny_ = opts.getInt("Ny",1);
    Nx_ = sites_.N()/Ny_;
    const Real J = opts.getReal("J",1.);
    Jz_ = opts.getReal("Jz",J);
    Jxy_ = opts.getReal("Jxy",J);
    boundary_h_ = opts.getReal("Boundary_h",0.);
    }

void inline TriangularHeisenberg::
init_()
    {
    if(initted_) return;

    H = MPO(sites_);

    const int Ns = sites_.N();
    const int nop = 3;
    const int max_mpo_dist = Ny_+1;
    const int k = 5+nop*(max_mpo_dist-1);

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;

    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Anc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(sites_.si(n),sites_.siP(n),row,col);

        W += sites_.op("Id",n) * row(1) * col(1);
        W += sites_.op("Id",n) * row(k) * col(k);

        W += sites_.op("Sz",n) * row(2) * col(1);
        W += sites_.op("S+",n) * row(3) * col(1);
        W += sites_.op("S-",n) * row(4) * col(1);

        //Horizontal bonds, connect n -> n+Ny_
        int mpo_dist = Ny_; 
        W += sites_.op("Sz",n) * row(k) * col(2+nop*(mpo_dist-1)) * Jz_;
        W += sites_.op("S-",n) * row(k) * col(3+nop*(mpo_dist-1)) * Jxy_/2;
        W += sites_.op("S+",n) * row(k) * col(4+nop*(mpo_dist-1)) * Jxy_/2;

        //Add boundary field if requested
        const int x = (n-1)/Ny_+1, y = (n-1)%Ny_+1;
        if(boundary_h_ != 0 && (x == 1 || x == Nx_))
            {
            Real eff_h = boundary_h_;
            eff_h *= ((x+y-2)%3==0 ? -1 : 0.5);
            printfln("Applying a pinning field of %.2f at site %d (%d,%d)\n",eff_h,n,x,y);
            W += sites_.op("Sz",n) * ITensor(row(k),col(1)) * eff_h;
            }

        //The following only apply if Ny_ > 1:

        //String of identity ops
        for(int q = 1; q <= nop*(max_mpo_dist-1); ++q)
            { W += sites_.op("Id",n) * row(1+nop+q) * col(1+q); }

        //Square lattice periodic BC bond
        //Connects n -> n+(Ny_-1)
        if(y == 1)
            {
            mpo_dist = Ny_-1; 
            W += sites_.op("Sz",n) * row(k) * col(2+nop*(mpo_dist-1)) * Jz_;
            W += sites_.op("S-",n) * row(k) * col(3+nop*(mpo_dist-1)) * Jxy_/2;
            W += sites_.op("S+",n) * row(k) * col(4+nop*(mpo_dist-1)) * Jxy_/2;
            }

        //N.N. bond along column
        W += sites_.op("Sz",n) * row(k) * col(2) * Jz_;
        W += sites_.op("S-",n) * row(k) * col(3) * Jxy_/2;
        W += sites_.op("S+",n) * row(k) * col(4) * Jxy_/2;

        //Diagonal bonds
        if(y != Ny_)
            {
            mpo_dist = Ny_+1; 
            W += sites_.op("Sz",n) * row(k) * col(2+nop*(mpo_dist-1)) * Jz_;
            W += sites_.op("S-",n) * row(k) * col(3+nop*(mpo_dist-1)) * Jxy_/2;
            W += sites_.op("S+",n) * row(k) * col(4+nop*(mpo_dist-1)) * Jxy_/2;
            }
        }

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

}; //namespace itensor


#endif
