//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_ISING_H
#define __ITENSOR_HAMS_ISING_H
#include "../mpo.h"
#include "../hams.h"

class Ising : public MPOBuilder
{
public:
    typedef MPOBuilder Parent;

    //
    // Constructors
    //

    Ising(const Model& model_);

    Ising(const Model& model_, Real hx, int ny = 1);

    //
    // Accessor Methods
    //

    Real
    J() const { return J_; }
    void
    J(Real val) { initted = false; J_ = val; }

    Real
    hx() const { return hx_; }
    void
    hx(Real val) { initted = false; hx_ = val; }

    int
    ny() const { return ny_; }
    void
    ny(int val) { initted = false; ny_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

private:

    const Model& model;
    int ny_,nx_;
    Real J_, hx_;
    bool initted;
    MPO H;

    void init_();

}; //class Ising

inline Ising::
Ising(const Model& model_) 
    : Parent(model_), 
      model(model_), 
      ny_(1),
      nx_(this->Ns),
      J_(1), 
      hx_(0),
      initted(false)
    { }

inline Ising::
Ising(const Model& model_, Real hx, int ny)
    : Parent(model_), 
      model(model_), 
      ny_(ny),
      nx_(this->Ns/ny_),
      J_(1), 
      hx_(hx),
      initted(false)
    { }

void inline Ising::
init_()
    {
    if(initted) return;

    H = MPO(model);

    const int max_mpo_dist = ny_;
    const int k = 3+(max_mpo_dist-1);

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;

    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.AAnc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model.si(n),model.siP(n),row,col);

        W += model.id(n) * row(1) * col(1);
        W += model.id(n) * row(k) * col(k);

        W += model.sz(n) * row(2) * col(1);

        //Transverse field
        if(hx_ != 0)
            W += model.sx(n) * ITensor(row(k),col(1)) * hx_;

        //Horizontal bonds (N.N in 1d)
        int mpo_dist = ny_; 
        W += model.sz(n) * row(k) * col(2+(mpo_dist-1)) * J_;

        //
        //The following only apply if ny_ > 1:
        //

        //String of identity ops
        for(int q = 1; q <= (max_mpo_dist-1); ++q)
            { W += model.id(n) * row(2+q) * col(1+q); }

        //Periodic BC bond
        const int y = (n-1)%ny_+1;
        if(y == 1 && ny_ > 2)
            {
            int mpo_dist = ny_-1; 
            W += model.sz(n) * row(k) * col(2+(mpo_dist-1)) * J_;
            }

        //N.N. bond along column
        if(y != ny_)
            {
            W += model.sz(n) * row(k) * col(2) * J_;
            }
        }

        H.AAnc(1) *= makeLedge(links.at(0));
        H.AAnc(Ns) *= makeRedge(links.at(Ns));

    initted = true;
    }


#endif
