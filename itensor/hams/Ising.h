//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_ISING_H
#define __ITENSOR_HAMS_ISING_H
#include "../mpo.h"
#include "../hams.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

class Ising
    {
    public:

    //
    // Constructors
    //

    Ising(const Model& model,
          const OptSet& opts = Global::opts());

    //
    // Accessor Methods
    //

    Real
    J() const { return J_; }
    void
    J(Real val) { initted_ = false; J_ = val; }

    Real
    hx() const { return hx_; }
    void
    hx(Real val) { initted_ = false; hx_ = val; }

    int
    Ny() const { return Ny_; }
    void
    Ny(int val) { initted_ = false; Ny_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

private:

    const Model& model_;
    int Ny_,
        Nx_;
    Real J_, 
         hx_;
    bool initted_;
    MPO H;

    void init_();

    }; //class Ising

inline Ising::
Ising(const Model& model,
      const OptSet& opts)
    :
    model_(model), 
    initted_(false)
    { 
    Ny_ = opts.getInt("Ny",1);
    Nx_ = model_.NN()/Ny_;
    J_  = opts.getReal("J",1.);
    hx_  = opts.getReal("hx",0.);
    }


void inline Ising::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.NN();
    const int max_mpo_dist = Ny_;
    const int k = 3+(max_mpo_dist-1);

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Aref(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        W += model_.id(n) * row(1) * col(1);
        W += model_.id(n) * row(k) * col(k);

        W += model_.sz(n) * row(2) * col(1);

        //Transverse field
        if(hx_ != 0)
            {
            W += model_.sx(n) * row(k) * col(1) * hx_;
            }

        //Horizontal bonds (N.N in 1d)
        int mpo_dist = Ny_; 
        W += model_.sz(n) * row(k) * col(2+(mpo_dist-1)) * J_;

        //
        //The following only apply if ny_ > 1:
        //

        //String of identity ops
        for(int q = 1; q <= (max_mpo_dist-1); ++q)
            { W += model_.id(n) * row(2+q) * col(1+q); }

        //Periodic BC bond
        const int y = (n-1)%Ny_+1;
        if(y == 1 && Ny_ > 2)
            {
            int mpo_dist = Ny_-1; 
            W += model_.sz(n) * row(k) * col(2+(mpo_dist-1)) * J_;
            }

        //N.N. bond along column
        if(y != Ny_)
            {
            W += model_.sz(n) * row(k) * col(2) * J_;
            }
        }

    H.Aref(1) *= ITensor(links.at(0)(k));
    H.Aref(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

#undef Cout
#undef Endl
#undef Format

#endif
