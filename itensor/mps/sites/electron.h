//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HUBBARD_H
#define __ITENSOR_HUBBARD_H
#include "itensor/mps/siteset.h"

namespace itensor {

class ElectronSite;

using Electron = BasicSiteSet<ElectronSite>;

class ElectronSite
    {
    Index s;
    public:

    ElectronSite() { }

    ElectronSite(Index I) : s(I) { }

    ElectronSite(int n, Args const& args = Args::global())
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveNf = args.getBool("ConserveNf",conserveQNs);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);
        int Up = (conserveSz ? +1 : 0),
            Dn = -Up;
        auto ts = format("Site,Elec,n=%d",n);
        if(conserveNf)
            {
            s = Index{QN({"Sz", 0},{"Nf",0,-1}),1,
                      QN({"Sz",Up},{"Nf",1,-1}),1,
                      QN({"Sz",Dn},{"Nf",1,-1}),1,
                      QN({"Sz", 0},{"Nf",2,-1}),1,Out,ts};
            }
        else //don't conserve Nf, only fermion parity
            {
            if(!conserveSz) Error("One of ConserveSz or ConserveNf must be true for Electron sites");
            s = Index{QN({"Sz", 0},{"Pf",0,-1}),1,
                      QN({"Sz",+1},{"Pf",1,-1}),1,
                      QN({"Sz",-1},{"Pf",1,-1}),1,
                      QN({"Sz", 0},{"Pf",0,-1}),1,Out,ts};
            }
        }

    Index
    index() const { return s; }

    IndexVal
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
        if(state == "S" || state == "UpDn") 
            {
            return s(4);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        IndexVal Em(s(1)),
                 EmP(sP(1)),
                 Up(s(2)),
                 UpP(sP(2)),
                 Dn(s(3)),
                 DnP(sP(3)),
                 UD(s(4)),
                 UDP(sP(4));

        ITensor Op(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(Up,UpP,1);
            Op.set(UD,UDP,1);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(Dn,DnP,1);
            Op.set(UD,UDP,1);
            }
        else
        if(opname == "Nupdn")
            {
            Op.set(UD,UDP,1);
            }
        else
        if(opname == "Ntot")
            {
            Op.set(Up,UpP,1);
            Op.set(Dn,DnP,1);
            Op.set(UD,UDP,2);
            }
        else
        if(opname == "Cup")
            {
            Op.set(Up,EmP,1); 
            Op.set(UD,DnP,1); 
            }
        else
        if(opname == "Cdagup")
            {
            Op.set(Em,UpP,1); 
            Op.set(Dn,UDP,1);
            }
        else
        if(opname == "Cdn")
            {
            Op.set(Dn,EmP,1); 
            Op.set(UD,UpP,-1); 
            }
        else
        if(opname == "Cdagdn")
            {
            Op.set(Em,DnP,1); 
            Op.set(Up,UDP,-1);
            }
        else
        if(opname == "Aup")
            {
            Op.set(Up,EmP,1); 
            Op.set(UD,DnP,1); 
            }
        else
        if(opname == "Adagup")
            {
            Op.set(Em,UpP,1); 
            Op.set(Dn,UDP,1);
            }
        else
        if(opname == "Adn")
            {
            Op.set(Dn,EmP,1); 
            Op.set(UD,UpP,1); 
            }
        else
        if(opname == "Adagdn")
            {
            Op.set(Em,DnP,1); 
            Op.set(Up,UDP,1);
            }
        else
        if(opname == "FermiPhase" || opname == "F")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,-1);
            Op.set(Dn,DnP,-1);
            Op.set(UD,UDP,+1);
            }
        else
        if(opname == "Fup")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,-1);
            Op.set(Dn,DnP,+1);
            Op.set(UD,UDP,-1);
            }
        else
        if(opname == "Fdn")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,+1);
            Op.set(Dn,DnP,-1);
            Op.set(UD,UDP,-1);
            }
        else
        if(opname == "Sz")
            {
            Op.set(Up,UpP,+0.5); 
            Op.set(Dn,DnP,-0.5);
            }
        else
        if(opname == "S+")
            {
            Op.set(Dn,UpP,1); 
            }
        else
        if(opname == "S-")
            {
            Op.set(Up,DnP,1); 
            }
        else
        if(opname == "S2")
            {
            //S dot S on-site
            Op.set(Up,UpP,0.75); 
            Op.set(Dn,DnP,0.75);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

//
// Deprecated, for backwards compatability
//

using HubbardSite = ElectronSite;

using Hubbard = BasicSiteSet<HubbardSite>;


} //namespace itensor

#endif
