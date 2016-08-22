//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINHALF_H
#define __ITENSOR_SPINHALF_H
#include "itensor/mps/siteset.h"

namespace itensor {

class SpinHalf : public SiteSet
    {
    public:

    SpinHalf() { }

    SpinHalf(int N);

    private:

    IQIndexVal
    getState(int i, const String& state) const;

    IQTensor
    getOp(int i, const String& opname, const Args& args) const;

    DefaultOpsT
    getDefaultOps(const Args& args) const;

    virtual void
    constructSites();

    //Data members -----------------

    static DefaultOpsT
    initDefaultOps()
        {
        DefaultOpsT dops;
        dops.push_back("Sz");
        return dops;
        }
        
    };

inline SpinHalf::
SpinHalf(int N)
    : SiteSet(N)
    { 
    constructSites();
    }

inline void SpinHalf::
constructSites()
    {
    for(int j = 1; j <= N(); ++j)
        {
        set(j,{nameint("S=1/2 ",j),
               Index(nameint("Up ",j),1,Site),QN("Sz=",+1),
               Index(nameint("Dn ",j),1,Site),QN("Sz=",-1)});
        }
    }

inline IQIndexVal SpinHalf::
getState(int i, const String& state) const
    {
    if(state == "Up") 
        {
        return si(i)(1);
        }
    else 
    if(state == "Dn") 
        {
        return si(i)(2);
        }
    else
        {
        Error("State " + state + " not recognized");
        }
    return si(i)(1);
    }

inline IQTensor SpinHalf::
getOp(int i, const String& opname, const Args& args) const
    {
    auto s  = si(i);
    auto sP = prime(s);

    auto Up = s(1);
    auto UpP = sP(1);
    auto Dn = s(2);
    auto DnP = sP(2);

    auto Op = IQTensor(dag(s),sP);

    if(opname == "Sz")
        {
        Op.set(Up,UpP,+0.5);
        Op.set(Dn,DnP,-0.5);
        }
    else
    if(opname == "Sx")
        {
        //mixedIQTensor call needed here
        //because as an IQTensor, Op would
        //not have a well defined QN flux
        Op = mixedIQTensor(s,sP);
        Op.set(Up,DnP,+0.5);
        Op.set(Dn,UpP,+0.5);
        }
    else
    if(opname == "ISy")
        {
        //mixedIQTensor call needed here
        //because as an IQTensor, Op would
        //not have a well defined QN flux
        Op = mixedIQTensor(s,sP);
        Op.set(Up,DnP,-0.5);
        Op.set(Dn,UpP,+0.5);
        }
    else
    if(opname == "Sy")
        {
        //mixedIQTensor call needed here
        //because as an IQTensor, Op would
        //not have a well defined QN flux
        Op = mixedIQTensor(s,sP);
        Op.set(Up,DnP,+0.5*Cplx_i);
        Op.set(Dn,UpP,-0.5*Cplx_i);
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
    if(opname == "projUp")
        {
        Op.set(Up,UpP,1);
        }
    else
    if(opname == "projDn")
        {
        Op.set(Dn,DnP,1);
        }
    else
    if(opname == "S2")
        {
        Op.set(Up,UpP,0.75);
        Op.set(Dn,DnP,0.75);
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

SpinHalf::DefaultOpsT inline SpinHalf::
getDefaultOps(Args const& args) const
    {
    static const auto dops_ = initDefaultOps();
    return dops_;
    }

} //namespace itensor

#endif
