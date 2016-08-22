//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HUBBARD_H
#define __ITENSOR_HUBBARD_H
#include "itensor/mps/siteset.h"

namespace itensor {

class Hubbard : public SiteSet
    {
    public:

    Hubbard();

    Hubbard(int N, 
            Args const& args = Global::args());

    bool
    conserveNf() const { return conserveNf_; }

    private:

    void
    doRead(std::istream& s);

    virtual IQIndexVal
    getState(int i, String const& state) const;

    virtual IQTensor
    getOp(int i, String const& opname, Args const& args = Global::args()) const;

    DefaultOpsT
    getDefaultOps(Args const& args) const;

    void 
    constructSites(int N);

        
    //Data members -----------------

    bool conserveNf_,
         conserveSz_;

    static DefaultOpsT
    initDefaultOps()
        {
        DefaultOpsT dops;
        dops.push_back("Ntot");
        dops.push_back("Nup");
        dops.push_back("Ndn");
        dops.push_back("Nupdn");
        dops.push_back("S2");
        return dops;
        }

    };

inline Hubbard::
Hubbard()
  : conserveNf_(true),
    conserveSz_(true)
    { }

inline Hubbard::
Hubbard(int N, Args const& args)
  : SiteSet(N)
    { 
    conserveNf_ = args.getBool("ConserveNf",true);
    conserveSz_ = args.getBool("ConserveSz",true);
    constructSites(N);
    }

void inline Hubbard::
constructSites(int N)
    {
    int Up = (conserveSz_ ? +1 : 0),
        Dn = -Up;
    if(conserveNf_)
        {
        for(auto j : range1(N))
            {
            set(j,{nameint("site=",j),
                Index(nameint("Emp ",j),1,Site), QN("Sz=", 0,"Nf=",0),
                Index(nameint("Up ",j),1,Site),  QN("Sz=",Up,"Nf=",1),
                Index(nameint("Dn ",j),1,Site),  QN("Sz=",Dn,"Nf=",1),
                Index(nameint("UpDn ",j),1,Site),QN("Sz=", 0,"Nf=",2)});
            }
        }
    else //don't conserve Nf, only fermion parity
        {
        if(!conserveSz_) Error("One of ConserveSz or ConserveNf must be true for Hubbard sites");

        for(auto j : range1(N))
            {
            set(j,{nameint("site=",j),
                Index(nameint("Emp ",j),1,Site), QN("Sz=", 0,"Pf=",0),
                Index(nameint("Up ",j),1,Site),  QN("Sz=",+1,"Pf=",1),
                Index(nameint("Dn ",j),1,Site),  QN("Sz=",-1,"Pf=",1),
                Index(nameint("UpDn ",j),1,Site),QN("Sz=", 0,"Pf=",0)});
            }
        }
    }

void inline Hubbard::
doRead(std::istream& s)
    {
    conserveNf_ = (si(1).qn(2)(2) == 1);
    }

inline IQIndexVal Hubbard::
getState(int i, String const& state) const
    {
    if(state == "0" || state == "Emp") 
        {
        return si(i)(1);
        }
    else 
    if(state == "+" || state == "Up") 
        {
        return si(i)(2);
        }
    else 
    if(state == "-" || state == "Dn") 
        {
        return si(i)(3);
        }
    else 
    if(state == "S" || state == "UpDn") 
        {
        return si(i)(4);
        }
    else
        {
        Error("State " + state + " not recognized");
        }
    return IQIndexVal{};
    }


inline IQTensor Hubbard::
getOp(int i, String const& opname, Args const& args) const
    {
    auto s = si(i);
    auto sP = prime(s);

    IQIndexVal Em(s(1)),
               EmP(sP(1)),
               Up(s(2)),
               UpP(sP(2)),
               Dn(s(3)),
               DnP(sP(3)),
               UD(s(4)),
               UDP(sP(4));

    IQTensor Op(dag(s),sP);

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
    if(opname == "S2")
        {
        //S dot S on-site
        Op.set(Up,UpP,0.75); 
        Op.set(Dn,DnP,0.75);
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

Hubbard::DefaultOpsT inline Hubbard::
getDefaultOps(Args const& args) const
    {
    static const auto dops_ = initDefaultOps();
    return dops_;
    }

} //namespace itensor

#endif
