#ifndef __ITENSOR_HAMS_TRIANGULARHEISENBERG_H
#define __ITENSOR_HAMS_TRIANGULARHEISENBERG_H
#include "../mpo.h"
#include "../hams.h"

class TriangularHeisenberg : public MPOBuilder
    {
    public:

    typedef MPOBuilder 
    Parent;

    TriangularHeisenberg(const Model& model_, int ny);

    TriangularHeisenberg(const Model& model_, int ny, Real boundary_h);

    Real
    J() const { return J_; }

    void
    J(Real val) { initted = false; J_ = val; }

    int
    ny() const { return ny_; }

    void
    ny(int val) { initted = false; ny_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    //------------------------------------------------------//

    private:

    const Model& model;
    int ny_,nx_;
    Real J_, boundary_h_;
    bool initted;
    MPO H;

    void init_();

    }; //class TriangularHeisenberg

inline TriangularHeisenberg::
TriangularHeisenberg(const Model& model_, int ny) 
    : Parent(model_), 
      model(model_), 
      ny_(ny),
      nx_(this->Ns/ny_),
      J_(1), 
      boundary_h_(0),
      initted(false)
    { }

inline TriangularHeisenberg::
TriangularHeisenberg(const Model& model_, int ny, Real boundary_h) 
    : Parent(model_), 
      model(model_), 
      ny_(ny),
      nx_(this->Ns/ny_),
      J_(1), 
      boundary_h_(boundary_h),
      initted(false)
    { }

void inline TriangularHeisenberg::
init_()
    {
    if(initted) return;

    H = MPO(model);

    const int nop = 3;
    const int max_mpo_dist = ny_+1;
    const int k = 5+nop*(max_mpo_dist-1);

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
        W += model.sp(n) * row(3) * col(1);
        W += model.sm(n) * row(4) * col(1);

        //Horizontal bonds, connect n -> n+ny_
        int mpo_dist = ny_; 
        W += model.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * J_;
        W += model.sm(n) * row(k) * col(3+nop*(mpo_dist-1)) * J_/2;
        W += model.sp(n) * row(k) * col(4+nop*(mpo_dist-1)) * J_/2;

        //Add boundary field if requested
        const int x = (n-1)/ny_+1, y = (n-1)%ny_+1;
        if(boundary_h_ != 0 && (x == 1 || x == nx_))
            {
            Real eff_h = boundary_h_;
            eff_h *= ((x+y-2)%3==0 ? -1 : 0.5);
            std::cout << boost::format("Applying a pinning field of %.2f at site %d (%d,%d)\n")%eff_h%n%x%y;
            W += model.sz(n) * ITensor(row(k),col(1)) * eff_h;
            }

        //The following only apply if ny_ > 1:

        //String of identity ops
        for(int q = 1; q <= nop*(max_mpo_dist-1); ++q)
            { W += model.id(n) * row(1+nop+q) * col(1+q); }

        //Square lattice periodic BC bond
        //Connects n -> n+(ny_-1)
        if(y == 1)
            {
            mpo_dist = ny_-1; 
            W += model.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * J_;
            W += model.sm(n) * row(k) * col(3+nop*(mpo_dist-1)) * J_/2;
            W += model.sp(n) * row(k) * col(4+nop*(mpo_dist-1)) * J_/2;
            }

        //N.N. bond along column
        W += model.sz(n) * row(k) * col(2) * J_;
        W += model.sm(n) * row(k) * col(3) * J_/2;
        W += model.sp(n) * row(k) * col(4) * J_/2;

        //Diagonal bonds
        if(y != ny_)
            {
            mpo_dist = ny_+1; 
            W += model.sz(n) * row(k) * col(2+nop*(mpo_dist-1)) * J_;
            W += model.sm(n) * row(k) * col(3+nop*(mpo_dist-1)) * J_/2;
            W += model.sp(n) * row(k) * col(4+nop*(mpo_dist-1)) * J_/2;
            }
        }

    H.AAnc(1) *= makeLedge(links.at(0));
    H.AAnc(Ns) *= makeRedge(links.at(Ns));

    initted = true;
    }


#endif
