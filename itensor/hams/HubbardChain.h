//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_HUBBARDCHAIN_H
#define __ITENSOR_HAMS_HUBBARDCHAIN_H
#include "../mpo.h"
#include "../model/hubbard.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

class HubbardChain
    {
    public:

    HubbardChain(const Hubbard& model,
                 const OptSet& opts = Global::opts());

    operator MPO() { init_(); return H.toMPO(); }

    operator IQMPO() { init_(); return H; }

    private:

    ///////////////////
    //
    // Data Members

    const Hubbard& model_;
    Real t_,U_;
    bool initted_;
    bool infinite_;
    IQMPO H;

    //
    //////////////////

    void 
    init_();

    }; //class HubbardChain

inline HubbardChain::
HubbardChain(const Hubbard& model, 
             const OptSet& opts)
    : 
    model_(model), 
    initted_(false)
    { 
    U_ = opts.getReal("U",0);
    t_ = opts.getReal("t",1);
    infinite_ = opts.getBool("Infinite",false);
    }

void inline HubbardChain::
init_()
    {
    if(initted_) return;

    H = IQMPO(model_);

    const int Ns = model_.N();

    std::vector<IQIndex> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) 
        {
        links.at(l) = IQIndex(nameint("Hl",l),
                              Index(nameint("h00_",l),2),QN( 0, 0),
                              Index(nameint("hup_",l),1),QN(+1,+1),
                              Index(nameint("hum_",l),1),QN(+1,-1),
                              Index(nameint("hdp_",l),1),QN(-1,+1),
                              Index(nameint("hdm_",l),1),QN(-1,-1)
                              );
        }

    const
    IQIndex last = (infinite_ ? links.at(0) : links.at(Ns));

    for(int n = 1; n <= Ns; ++n)
        {
        IQTensor& W = H.Anc(n);
        const
        IQIndex row = conj(links[n-1]), 
                col = (n==Ns ? last : links[n]);

        W = IQTensor(conj(model_.si(n)),model_.siP(n),row,col);

        //Identity strings
        W += model_.op("Id",n) * row(1) * col(1);
        W += model_.op("Id",n) * row(2) * col(2);

        //Hubbard U
        W += model_.op("Nupdn",n) * row(2) * col(1) * U_;

        //Kinetic energy/hopping terms, defined as -t_*(c^d_i c_{i+1} + h.c.)
        W += model_.op("Aup*F",n)    *row(2)*col(3)*(+t_);
        W += model_.op("Adagdn",n)   *row(2)*col(4)*(-t_);
        W += model_.op("Adn",n)      *row(2)*col(5)*(+t_);
        W += model_.op("Adagup*F",n) *row(2)*col(6)*(-t_);

        W += model_.op("Adagup",n)  *row(3)*col(1);
        W += model_.op("F*Adn",n)   *row(4)*col(1);
        W += model_.op("F*Adagdn",n)*row(5)*col(1);
        W += model_.op("Aup",n)     *row(6)*col(1);
        }

    const
    IQTensor LH(links.at(0)(2)),
             RH(conj(last)(1)); 

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

#undef Cout
#undef Endl
#undef Format

#endif
