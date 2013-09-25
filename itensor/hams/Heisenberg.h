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

    operator MPO() { init_(); return H.toMPO(); }

    operator IQMPO() { init_(); return H; }

    private:

    //////////////////
    //
    // Data Members

    const Model& model_;
    int Ny_,
        Nx_;
    Real J_, 
         Boundary_h_,
         Jz_;
    bool initted_;
    IQMPO H;

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
    Jz_ = opts.getReal("Jz",J_);
    }


void inline Heisenberg::
init_()
    {
    if(initted_) return;

    H = IQMPO(model_);

    const int Ns = model_.N(),
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

    for(int n = 1; n <= Ns; ++n)
        {
        IQTensor& W = H.Anc(n);
        IQIndex &row = links.at(n-1), 
                &col = links.at(n);

        W = IQTensor(conj(model_.si(n)),model_.siP(n),conj(row),col);

        W += model_.op("Id",n) * row(1) * col(1); //ending state
        W += model_.op("Id",n) * row(2) * col(2); //starting state

        //
        //Sz Sz terms
        //
        W += model_.op("Sz",n) * row(3) * col(1);
        //Horizontal bonds
        W += model_.op("Sz",n) * row(start) * col(kd) * Jz_;

        //
        //S+ S- terms
        //
        W += model_.op("Sm",n) * row(kd+1) * col(1);
        //Horizontal bonds
        W += model_.op("Sp",n) * row(start) * col(kd+kpm) * J_/2;

        //
        //S- S+ terms
        //
        W += model_.op("Sp",n) * row(kd+kpm+1) * col(1);
        //Horizontal bonds
        W += model_.op("Sm",n) * row(start) * col(kd+2*kpm) * J_/2;

        //Add boundary field if requested
        const int x = (n-1)/Ny_+1, y = (n-1)%Ny_+1;
        if(Boundary_h_ != 0 && (x == 1 || x == Nx_))
            {
            Real eff_h = Boundary_h_;
            if(J_ > 0) eff_h *= (x%2==1 ? -1 : 1)*(y%2==1 ? -1 : 1);
            //cerr << format("Doing a staggered bf of %.2f at site %d (%d,%d)\n")%eff_h%n%x%y;
            W += model_.op("Sz",n) * row(start) * col(1) * eff_h;
            }

        //The following only apply if Ny_ > 1:

        //Strings of off-diagonal identity ops
        for(int q = 1; q <= Ny_-1; ++q)
            { 
            W += model_.op("Id",n) * row(3+q) * col(2+q); 
            W += model_.op("Id",n) * row(kd+1+q) * col(kd+q); 
            W += model_.op("Id",n) * row(kd+kpm+1+q) * col(kd+kpm+q); 
            }

        //Periodic BC bond (only for width 3 ladders or greater)
        if(y == 1 && Ny_ >= 3)
            {
            W += model_.op("Sz",n) * row(start) * col(kd-1) * Jz_;
            W += model_.op("Sp",n) * row(start) * col(kd+kpm-1) * J_/2;
            W += model_.op("Sm",n) * row(start) * col(kd+2*kpm-1) * J_/2;
            }

        //N.N. bond along column
        if(y != Ny_)
            {
            W += model_.op("Sz",n) * row(start) * col(3) * J_;
            W += model_.op("Sp",n) * row(start) * col(kd+1) * J_/2;
            W += model_.op("Sm",n) * row(start) * col(kd+kpm+1) * J_/2;
            }
        }

    H.Anc(1) *= IQTensor(links.at(0)(start));
    H.Anc(Ns) *= IQTensor(conj(links.at(Ns))(1));

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
        Index &row = links.at(n-1), &col = links.at(n);

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        W += model_.id(n) * row(1) * col(1);
        W += model_.id(n) * row(k) * col(k);

        W += model_.sz(n) * row(2) * col(1);
        W += model_.sp(n) * row(3) * col(1);
        W += model_.sm(n) * row(4) * col(1);

        //Horizontal bonds
        int mpo_dist = Ny_; 
        W += model_.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * Jz_;
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

        //Periodic BC bond (only for width 3 ladders or greater)
        if(y == 1 && Ny_ >= 3)
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
*/

#undef Cout
#undef Endl
#undef Format

#endif
