#ifndef __ITENSOR_HAMS_HUBBARDCHAIN_H
#define __ITENSOR_HAMS_HUBBARDCHAIN_H
#include "../mpo.h"
#include "../hams.h"

class HubbardChain : public MPOBuilder
{
public:
    typedef MPOBuilder Parent;

    HubbardChain(const Hubbard::Model& model_);

    HubbardChain(const Hubbard::Model& model_, Real U);

    Real
    t() const { return t_; }

    void
    t(Real val) { initted = false; t_ = val; }

    Real
    U() const { return U_; }

    void
    U(Real val) { initted = false; U_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

private:

    const Hubbard::Model& model;
    Real t_,U_;
    bool initted;
    MPO H;

    void init_();

}; //class HubbardChain

inline HubbardChain::
HubbardChain(const Hubbard::Model& model_) 
    : Parent(model_), 
      model(model_), 
      t_(1), 
      U_(0), 
      initted(false)
    { }

inline HubbardChain::
HubbardChain(const Hubbard::Model& model_, Real U) 
    : Parent(model_), 
      model(model_), 
      t_(1), 
      U_(U), 
      initted(false)
    { }

void inline HubbardChain::
init_()
    {
    if(initted) return;

    H = MPO(model);

    const int k = 6;

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;
    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.AAnc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model.si(n),model.siP(n),row,col);

        W += model.NupNdn(n) * row(k) * col(1) * U_;
        W += multSiteOps(model.FermiPhase(n),model.Cup(n)) * row(k) * col(2) * t_;
        W += multSiteOps(model.FermiPhase(n),model.Cdn(n)) * row(k) * col(3) * t_;
        W += multSiteOps(model.Cdagup(n),model.FermiPhase(n)) * row(k) * col(4) * t_;
        W += multSiteOps(model.Cdagdn(n),model.FermiPhase(n)) * row(k) * col(5) * t_;
        W += model.id(n) * row(k) * col(k);

        W += model.id(n) * row(1) * col(1);
        W += model.Cdagup(n) * row(2) * col(1) * (-1.0);
        W += model.Cdagdn(n) * row(3) * col(1) * (-1.0);
        W += model.Cup(n) * row(4) * col(1) * (-1.0);
        W += model.Cdn(n) * row(5) * col(1) * (-1.0);
        }

    H.AAnc(1) *= makeLedge(links.at(0));
    H.AAnc(Ns) *= makeRedge(links.at(Ns));

    initted = true;
    }

#endif
