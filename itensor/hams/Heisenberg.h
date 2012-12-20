//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_HEISENBERG_H
#define __ITENSOR_HAMS_HEISENBERG_H

#include "../mpo.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

class Heisenberg
    {
    public:

    Heisenberg(const Model& model, 
               const OptSet& opts = Global::opts());

    Real
    J() const { return J_; }
    void
    J(Real val) { initted_ = false; J_ = val; }

    int
    Ny() const { return Ny_; }
    void
    Ny(int val) { initted_ = false; Ny_ = val; }

    Real
    Boundary_h() const { return Boundary_h_; }
    void
    Boundary_h(Real val) { initted_ = false; Boundary_h_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    //////////////////
    //
    // Data Members

    const Model& model_;
    int Ny_,
        Nx_;
    Real J_, 
         Boundary_h_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class Heisenberg

inline Heisenberg::
Heisenberg(const Model& model, const OptSet& opts)
    : 
    model_(model), 
    initted_(false)
    { 
    Ny_ = opts.getInt("Ny",1);
    Nx_ = model_.N()/Ny_;
    J_ = opts.getReal("J",1.);
    Boundary_h_ = opts.getReal("Boundary_h",0.);
    }

void inline Heisenberg::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.N();
    const int nop = 3;
    const int max_mpo_dist = Ny_;
    const int k = 5+nop*(max_mpo_dist-1);

    std::vector<Index> links(model_.N()+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Anc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        W += model_.id(n) * row(1) * col(1);
        W += model_.id(n) * row(k) * col(k);

        W += model_.sz(n) * row(2) * col(1);
        W += model_.sp(n) * row(3) * col(1);
        W += model_.sm(n) * row(4) * col(1);

        //Horizontal bonds
        int mpo_dist = Ny_; 
        W += model_.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * J_;
        W += model_.sm(n) * row(k) * col(3+nop*(mpo_dist-1)) * J_/2;
        W += model_.sp(n) * row(k) * col(4+nop*(mpo_dist-1)) * J_/2;

        //Add boundary field if requested
        const int x = (n-1)/Ny_+1, y = (n-1)%Ny_+1;
        if(Boundary_h_ != 0 && (x == 1 || x == Nx_))
            {
            Real eff_h = Boundary_h_;
            if(J_ > 0) eff_h *= (x%2==1 ? -1 : 1)*(y%2==1 ? -1 : 1);
            //cerr << format("Doing a staggered bf of %.2f at site %d (%d,%d)\n")%eff_h%n%x%y;
            W += model_.sz(n) * ITensor(row(k),col(1)) * eff_h;
            }

        //The following only apply if Ny_ > 1:

        //String of identity ops
        for(int q = 1; q <= nop*(max_mpo_dist-1); ++q)
            { W += model_.id(n) * row(1+nop+q) * col(1+q); }

        //Periodic BC bond
        if(y == 1 && Ny_ > 1)
            {
            int mpo_dist = Ny_-1; 
            W += model_.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * J_;
            W += model_.sm(n) * row(k) * col(3+nop*(mpo_dist-1)) * J_/2;
            W += model_.sp(n) * row(k) * col(4+nop*(mpo_dist-1)) * J_/2;
            }

        //N.N. bond along column
        if(y != Ny_)
            {
            W += model_.sz(n) * row(k) * col(2) * J_;
            W += model_.sm(n) * row(k) * col(3) * J_/2;
            W += model_.sp(n) * row(k) * col(4) * J_/2;
            }
        }

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

#undef Cout
#undef Endl
#undef Format

#endif
