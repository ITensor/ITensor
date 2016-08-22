//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINONE_H
#define __ITENSOR_SPINONE_H
#include "itensor/mps/siteset.h"

namespace itensor {

class SpinOne : public SiteSet
    {
    public:

    SpinOne() { }

    SpinOne(int N, 
            Args const& args = Args::global());

    private:

    void
    constructSites(Args const& args);

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const Args& args) const;

    };

inline SpinOne::
SpinOne(int N, const Args& args)
    : SiteSet(N)
    { 
    constructSites(args);
    }

inline void SpinOne::
constructSites(Args const& args)
    {
    for(int j = 1; j <= N(); ++j)
        {
        if((args.getBool("SHalfEdge",false) && (j==1 || j == N()))
           || (args.getBool("SHalfLeftEdge",false) && j==1))
            {
            if(args.getBool("Verbose",false))
                {
                println("Placing a S=1/2 at site ",j);
                }

            set(j,{nameint("S=1/2 site=",j),
                Index(nameint("Up:site",j),1,Site),QN("Sz=",+1),
                Index(nameint("Dn:site",j),1,Site),QN("Sz=",-1)});
            }
        else
            {
            set(j,{nameint("S=1 site=",j),
                Index(nameint("Up:site",j),1,Site),QN("Sz=",+2),
                Index(nameint("Z0:site",j),1,Site),QN("Sz=", 0),
                Index(nameint("Dn:site",j),1,Site),QN("Sz=",-2)});
            }
        }
    }

inline IQIndexVal SpinOne::
getState(int i, const String& state) const
    {
    int st = -1;
    if(state == "Up" || state == "+") 
        {
        st = 1;
        }
    else
    if(state == "Z0" || state == "0")
        {
        if(si(i).m() == 2) Error("Z0 not defined for spin 1/2 site");
        st = 2;
        }
    else
    if(state == "Dn" || state == "-")
        {
        st = si(i).m();
        }
    else
        {
        Error("State " + state + " not recognized");
        }
    return si(i)(st);
    }

inline IQTensor SpinOne::
getOp(int i, const String& opname, const Args& args) const
    {
    auto s = si(i);
    auto sP = prime(s);

    auto Up  = s(1);
    auto UpP = sP(1);
    auto Z0  = s(2);
    auto Z0P = sP(2);
    auto Dn  = s(s.m());
    auto DnP = sP(s.m());

    auto Op = IQTensor(dag(s),sP);

    if(opname == "Sz")
        {
        if(s.m() == 2)
            {
            Op.set(Up,UpP,+0.5);
            Op.set(Dn,DnP,-0.5);
            }
        else
            {
            Op.set(Up,UpP,+1.0);
            Op.set(Dn,DnP,-1.0);
            }
        }
    else
    if(opname == "Sx")
        {
        //mixedIQTensor call needed here
        //because as an IQTensor, Op would
        //not have a well defined QN flux
        Op = mixedIQTensor(s,sP);
        if(s.m() == 2)
            {
            Op.set(Up,DnP,+0.5);
            Op.set(Dn,UpP,+0.5);
            }
        else
            {
            Op.set(Up,Z0P,ISqrt2); 
            Op.set(Z0,UpP,ISqrt2);
            Op.set(Z0,DnP,ISqrt2); 
            Op.set(Dn,Z0P,ISqrt2);
            }
        }
    else
    if(opname == "ISy")
        {
        //mixedIQTensor call needed here
        //because as an IQTensor, Op would
        //not have a well defined QN flux
        Op = mixedIQTensor(s,sP);
        if(s.m() == 2)
            {
            Op.set(Up,DnP,-0.5);
            Op.set(Dn,UpP,+0.5);
            }
        else
            {
            Op.set(Up,Z0P,+ISqrt2); 
            Op.set(Z0,UpP,-ISqrt2);
            Op.set(Z0,DnP,+ISqrt2); 
            Op.set(Dn,Z0P,-ISqrt2);
            }
        }
    else
    if(opname == "Sy")
        {
        //mixedIQTensor call needed here
        //because as an IQTensor, Op would
        //not have a well defined QN flux
        Op = mixedIQTensor(s,sP);
        if(s.m() == 2)
            {
            Op.set(Up,DnP,+0.5_i);
            Op.set(Dn,UpP,-0.5_i);
            }
        else
            {
            Op.set(Up,Z0P,-ISqrt2*1_i); 
            Op.set(Z0,UpP,+ISqrt2*1_i);
            Op.set(Z0,DnP,-ISqrt2*1_i); 
            Op.set(Dn,Z0P,+ISqrt2*1_i);
            }
        }
    else
    if(opname == "Sp" || opname == "S+")
        {
        if(s.m() == 2)
            {
            Op.set(Dn,UpP,1);
            }
        else
            {
            Op.set(Dn,Z0P,Sqrt2);  
            Op.set(Z0,UpP,Sqrt2);
            }
        }
    else
    if(opname == "Sm" || opname == "S-")
        {
        if(s.m() == 2)
            {
            Op.set(Up,DnP,1);
            }
        else
            {
            Op.set(Up,Z0P,Sqrt2);
            Op.set(Z0,DnP,Sqrt2);
            }
        }
    else
    if(opname == "Sz2")
        {
        if(s.m() == 2) Error("Sz^2 only non-trivial for S=1 sites");
        Op.set(Up,UpP,1); 
        Op.set(Dn,DnP,1);
        }
    else
    if(opname == "Sx2")
        {
        if(s.m() == 2) Error("Sx^2 only non-trivial for S=1 sites");
        Op.set(Up,UpP,0.5); 
        Op.set(Up,DnP,0.5);
        Op.set(Z0,Z0P,1.0);
        Op.set(Dn,DnP,0.5); 
        Op.set(Dn,UpP,0.5);
        }
    else
    if(opname == "Sy2")
        {
        if(s.m() == 2) Error("Sy^2 only non-trivial for S=1 sites");
        Op.set(Up,UpP,+0.5); 
        Op.set(Up,DnP,-0.5);
        Op.set(Z0,Z0P,1);
        Op.set(Dn,DnP,+0.5); 
        Op.set(Dn,UpP,-0.5);
        }
    else
    if(opname == "projUp")
        {
        Op.set(Up,UpP,1);
        }
    else
    if(opname == "projZ0")
        {
        if(s.m() == 2) Error("Can only form projZ0 for S=1 sites");
        Op.set(Z0,Z0P,1);
        }
    else
    if(opname == "projDn")
        {
        Op.set(Dn,DnP,1);
        }
    else
    if(opname == "XUp")
        {
        //m = +1 state along x axis
        Op = IQTensor(s);
        Op.set(Up,0.5);
        Op.set(Z0,ISqrt2);
        Op.set(Dn,0.5);
        }
    else
    if(opname == "XZ0")
        {
        //m = 0 state along x axis
        Op = IQTensor(s);
        Op.set(Up,+ISqrt2);
        Op.set(Dn,-ISqrt2);
        }
    else
    if(opname == "XDn")
        {
        //m = -1 state along x axis
        Op = IQTensor(s);
        Op.set(Up,0.5);
        Op.set(Z0,-ISqrt2);
        Op.set(Dn,0.5);
        }
    else
    if(opname == "S2")
        {
        auto ssp1 = (s.m()==2 ? 0.75 : 2.);
        Op.set(Up,UpP,ssp1); 
        Op.set(Dn,DnP,ssp1);
        if(s.m() > 2)
            Op.set(Z0,Z0P,ssp1);
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

} //namespace itensor

#endif
