//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TJ_H
#define __ITENSOR_TJ_H
#include "itensor/mps/siteset.h"

namespace itensor {

class tJSite;

using tJ  = BasicSiteSet<tJSite>;

class tJSite
    {
    IQIndex s;
    public:

    tJSite() { }

    tJSite(IQIndex I) : s(I) { }

    tJSite(int n, Args const& args = Args::global())
        {
        s = IQIndex{nameint("tJ site=",n),
            Index(nameint("Emp ",n),1,Site), QN("Sz=", 0,"Nf=",0),
            Index(nameint("Up ",n),1,Site),  QN("Sz=",+1,"Nf=",1),
            Index(nameint("Dn ",n),1,Site),  QN("Sz=",-1,"Nf=",1)};
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "0" || state == "Emp") 
            {
            return s(1);
            }
        else 
        if(state == "+" || state == "Up") 
            {
            return s(2);
            }
        else 
        if(state == "-" || state == "Dn") 
            {
            return s(3);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal();
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        IQIndexVal Em(s(1)),
                   EmP(sP(1)),
                   Up(s(2)),
                   UpP(sP(2)),
                   Dn(s(3)),
                   DnP(sP(3));

        auto Op = IQTensor(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(Up,UpP,1);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(Dn,DnP,1);
            }
        else
        if(opname == "Ntot")
            {
            Op.set(Up,UpP,1);
            Op.set(Dn,DnP,1);
            }
        else
        if(opname == "Cup")
            {
            Op.set(Up,EmP,1); 
            }
        else
        if(opname == "Cdagup")
            {
            Op.set(Em,UpP,1); 
            }
        else
        if(opname == "Cdn")
            {
            Op.set(Dn,EmP,1); 
            }
        else
        if(opname == "Cdagdn")
            {
            Op.set(Em,DnP,1); 
            }
        else
        if(opname == "Aup")
            {
            Op.set(Up,EmP,1); 
            }
        else
        if(opname == "Adagup")
            {
            Op.set(Em,UpP,1); 
            }
        else
        if(opname == "Adn")
            {
            Op.set(Dn,EmP,1); 
            }
        else
        if(opname == "Adagdn")
            {
            Op.set(Em,DnP,1); 
            }
        else
        if(opname == "FermiPhase" || opname == "F")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,-1);
            Op.set(Dn,DnP,-1);
            }
        else
        if(opname == "Fup")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,-1);
            Op.set(Dn,DnP,+1);
            }
        else
        if(opname == "Fdn")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,+1);
            Op.set(Dn,DnP,-1);
            }
        else
        if(opname == "Sz")
            {
            Op.set(Up,UpP,+0.5); 
            Op.set(Dn,DnP,-0.5);
            }
        else
        if(opname == "Sx")
            {
            Op.set(Up,DnP,1); 
            Op.set(Dn,UpP,1);
            }
        else
        if(opname == "Sp" || opname == "S+")
            {
            Op.set(Dn,UpP,1);
            }
        else
        if(opname == "Sm" || opname == "S-")
            {
            Op.set(Up,DnP,1);
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
