//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_ISING_H
#define __ITENSOR_HAMS_ISING_H
#include "../mpo.h"

namespace itensor {

class Ising
    {
    public:

    //
    // Constructors
    //

    Ising(const SiteSet& sites,
          const OptSet& opts = Global::opts());

    operator MPO() { init(); return H; }

    private:

    ////////
    const SiteSet& sites_;
    int Ny_,
        Nx_;
    Real J_, 
         hx_;
    bool initted_;
    MPO H;
    bool infinite_;

    std::vector<Index> links;
    ////////

    void 
    init();

    }; //class Ising

inline Ising::
Ising(const SiteSet& sites,
      const OptSet& opts)
    :
    sites_(sites), 
    initted_(false)
    { 
    Ny_ = opts.getInt("Ny",1);
    Nx_ = sites_.N()/Ny_;
    J_  = opts.getReal("J",1.);
    hx_  = opts.getReal("hx",0.);
    infinite_ = opts.getBool("Infinite",false);
    }


void inline Ising::
init()
    {
    if(initted_) return;

    H = MPO(sites_);

    const int Ns = sites_.N();
    const int max_mpo_dist = Ny_;
    const int k = 3+(max_mpo_dist-1);

    links = std::vector<Index>(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    Index last = (infinite_ ? links.at(0) : links.at(Ns));

    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Anc(n);
        Index &row = links[n-1], 
              &col = (n==Ns ? last : links[n]);

        W = ITensor(sites_.si(n),prime(sites_.si(n)),row,col);

        W += sites_.op("Id",n) * row(1) * col(1);
        W += sites_.op("Id",n) * row(k) * col(k);

        W += sites_.op("Sz",n) * row(2) * col(1);

        //Transverse field
        if(hx_ != 0)
            {
            W += sites_.op("Sx",n) * row(k) * col(1) * (-hx_);
            }

        //Horizontal bonds (N.N in 1d)
        int mpo_dist = Ny_; 
        W += sites_.op("Sz",n) * row(k) * col(2+(mpo_dist-1)) * J_;

        //
        //The following only apply if Ny_ > 1:
        //

        //String of identity ops
        for(int q = 1; q <= (max_mpo_dist-1); ++q)
            { W += sites_.op("Id",n) * row(2+q) * col(1+q); }

        //Periodic BC bond
        const int y = (n-1)%Ny_+1;
        if(y == 1 && Ny_ > 2)
            {
            int mpo_dist = Ny_-1; 
            W += sites_.op("Sz",n) * row(k) * col(2+(mpo_dist-1)) * J_;
            }

        //N.N. bond along column
        if(y != Ny_)
            {
            W += sites_.op("Sz",n) * row(k) * col(2) * J_;
            }
        }

    const ITensor LH(links[0](k)),
                  RH(last(1));
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
}; //namespace itensor

#endif
