//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_Z3_H
#define __ITENSOR_Z3_H
#include "itensor/mps/siteset.h"

namespace itensor {

class Z3Site;

using Z3 = BasicSiteSet<Z3Site>;

class Z3Site
    {
    IQIndex s;
    public:

    Z3Site() { }

    Z3Site(IQIndex I) : s(I) { }

    Z3Site(int n, Args const& args = Args::global())
        {
        s = IQIndex{nameint("Z3 site=",n),
        Index(nameint("0|site",n),1,Site),QN({0,3}),
        Index(nameint("1|site",n),1,Site),QN({1,3}),
        Index(nameint("2|site",n),1,Site),QN({2,3})};
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "0") { return s(1); }
        else
        if(state == "1") { return s(2); }
        else
        if(state == "2") { return s(3); }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Zer = s(1);
        auto ZerP = sP(1);
        auto One = s(2);
        auto OneP = sP(2);
        auto Two = s(3);
        auto TwoP = sP(3);

        auto Op = IQTensor(dag(s),sP);

        if(opname == "N")
            {
            Op.set(One,OneP,1);
            Op.set(Two,TwoP,2);
            }
        else
        if(opname == "Sig")
            {
            Op.set(Zer,TwoP,1);
            Op.set(One,ZerP,1);
            Op.set(Two,OneP,1);
            }
        else
        if(opname == "SigDag")
            {
            Op.set(Two,ZerP,1);
            Op.set(Zer,OneP,1);
            Op.set(One,TwoP,1);
            }
        else
        if(opname == "Tau")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/3.)+sin(2.*Pi/3.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/3.)+sin(4.*Pi/3.)*1_i);
            }
        else
        if(opname == "TauDag")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/3.)-sin(2.*Pi/3.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/3.)-sin(4.*Pi/3.)*1_i);
            }
        else
        if(opname == "Proj0")
            {
            Op.set(Zer,ZerP,1);
            }
        else
        if(opname == "Proj1")
            {
            Op.set(One,OneP,1);
            }
        else
        if(opname == "Proj2")
            {
            Op.set(Two,TwoP,1);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

} //namespace itensor

#endif
