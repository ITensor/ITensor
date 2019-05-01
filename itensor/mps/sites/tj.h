//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TJ_H
#define __ITENSOR_TJ_H
#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class tJSite;

using tJ  = BasicSiteSet<tJSite>;

class tJSite
    {
    Index s;
    public:

    tJSite(Index I) : s(I) { }

    tJSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,tJ");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        if(args.getBool("ConserveQNs",true))
          {
          s = Index(QN({"Sz", 0},{"Nf",0}),1,
                    QN({"Sz",+1},{"Nf",1}),1,
                    QN({"Sz",-1},{"Nf",1}),1,Out,ts);
          }
        else
          {
          s = Index(3,ts);
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
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal();
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
                   DnP(sP(3));

        auto Op = ITensor(dag(s),sP);

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

    //
    // Deprecated, for backwards compatibility
    //

    tJSite(int n, Args const& args = Args::global())
        {
        *this = tJSite({args,"SiteNumber=",n});
        }

    };

} //namespace itensor

#endif
