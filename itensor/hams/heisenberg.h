#ifndef __ITENSOR_HAMS_HEISENBERG_H
#define __ITENSOR_HAMS_HEISENBERG_H
#include "../mpo.h"
#include "../hams.h"

class Heisenberg : public MPOBuilder
{
public:
    typedef MPOBuilder Parent;

    Heisenberg(const SpinModel& model_);

    Heisenberg(const SpinModel& model_, int ny);

    Heisenberg(const SpinModel& model_, int ny, Real boundary_h);

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

private:

    const SpinModel& model;
    int ny_,nx_;
    Real J_, boundary_h_;
    bool initted;
    MPO H;

    void init_();

}; //class Heisenberg

inline Heisenberg::
Heisenberg(const SpinModel& model_) 
    : Parent(model_), 
      model(model_), 
      ny_(1),
      nx_(this->Ns),
      J_(1), 
      boundary_h_(0),
      initted(false)
    { }

inline Heisenberg::
Heisenberg(const SpinModel& model_, int ny) 
    : Parent(model_), 
      model(model_), 
      ny_(ny),
      nx_(this->Ns/ny_),
      J_(1), 
      boundary_h_(0),
      initted(false)
    { }

inline Heisenberg::
Heisenberg(const SpinModel& model_, int ny, Real boundary_h) 
    : Parent(model_), 
      model(model_), 
      ny_(ny),
      nx_(this->Ns/ny_),
      J_(1), 
      boundary_h_(boundary_h),
      initted(false)
    { }

void inline Heisenberg::
init_()
    {
    if(initted) return;

    H = MPO(model);

    const int nop = 3;
    const int k = 2+nop*ny_;

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;

    /*
    if(model.periodic())
        {
        if(ny_ > 1) 
            Error("Periodic case not yet supported for ny > 1");
        if(boundary_h_ != 0)
            Error("Boundary field not allowed for periodic MPO");
        }
    //Make n == 1 tensor even for periodic, will replace below
    */

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

        W += model.sz(n) * row(k) * col(2+nop*(ny_-1)) * J_;
        W += model.sm(n) * row(k) * col(3+nop*(ny_-1)) * J_/2;
        W += model.sp(n) * row(k) * col(4+nop*(ny_-1)) * J_/2;

        //Add boundary field if requested

        const int x = (n-1)/ny_+1, y = (n-1)%ny_+1;
        if(boundary_h_ != 0 && (x == 1 || x == nx_))
            {
            Real eff_h = boundary_h_;
            if(J_ > 0) eff_h *= (x%2==1 ? -1 : 1)*(y%2==1 ? -1 : 1);
            //cerr << format("Doing a staggered bf of %.2f at site %d (%d,%d)\n")%eff_h%n%x%y;
            W += model.sz(n) * ITensor(row(k),col(1)) * eff_h;
            }

        //The following only applies if ny_ > 1

        for(int q = 1; q <= (ny_-1)*nop; ++q)
            { W += model.id(n) * row(1+nop+q) * col(1+q); }

        if(y == 1 && ny_ > 1)
            {
            W += model.sz(n) * row(k) * col(2+nop*(ny_-2)) * J_;
            W += model.sm(n) * row(k) * col(3+nop*(ny_-2)) * J_/2;
            W += model.sp(n) * row(k) * col(4+nop*(ny_-2)) * J_/2;
            }

        if(y != ny_)
            {
            W += model.sz(n) * row(k) * col(2) * J_;
            W += model.sm(n) * row(k) * col(3) * J_/2;
            W += model.sp(n) * row(k) * col(4) * J_/2;
            }
        }

    /*
    if(model.periodic())
        {
        ITensor& W = H.AAnc(1);
        Index &row = links[Ns], &col = links[1];

        W = ITensor(model.si(1),model.siP(1),row,col);

        W += model.id(1) * row(1) * col(k);

        W += model.sz(1) * row(2) * col(k);
        W += model.sp(1) * row(3) * col(k);
        W += model.sm(1) * row(4) * col(k);

        W += model.sz(1) * row(1) * col(2) * J_;
        W += model.sm(1) * row(1) * col(3) * J_/2;
        W += model.sp(1) * row(1) * col(4) * J_/2;
        }
    else //open BCs
        {
        */
        H.AAnc(1) *= makeLedge(links.at(0));
        H.AAnc(Ns) *= makeRedge(links.at(Ns));
        //}

    initted = true;
    }


#endif
