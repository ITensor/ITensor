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

    SpinHalf();

    SpinHalf(int N);

    private:

    int
    getN() const;

    const IQIndex&
    getSi(int i) const;

    IQIndexVal
    getState(int i, const String& state) const;

    IQTensor
    getOp(int i, const String& opname, const Args& args) const;

    DefaultOpsT
    getDefaultOps(const Args& args) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();

    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;

    static DefaultOpsT
    initDefaultOps()
        {
        DefaultOpsT dops;
        dops.push_back("Sz");
        return dops;
        }
        
    };

inline SpinHalf::
SpinHalf()
    : N_(-1)
    { }

inline SpinHalf::
SpinHalf(int N)
    : N_(N),
      site_(N_+1)
    { 
    constructSites();
    }

inline void SpinHalf::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("S=1/2 ",j),
            Index(nameint("Up ",j),1,Site),spin(+1),
            Index(nameint("Dn ",j),1,Site),spin(-1));
        }
    }

inline void SpinHalf::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

inline void SpinHalf::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

inline int SpinHalf::
getN() const
    { return N_; }

inline const IQIndex& SpinHalf::
getSi(int i) const
    { return site_.at(i); }

inline IQIndexVal SpinHalf::
getState(int i, const String& state) const
    {
    if(state == "Up") 
        {
        return getSi(i)(1);
        }
    else 
    if(state == "Dn") 
        {
        return getSi(i)(2);
        }
    else
        {
        Error("State " + state + " not recognized");
        return getSi(i)(1);
        }
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
getDefaultOps(const Args& args) const
    {
    static const std::vector<String> dops_(initDefaultOps());
    return dops_;
    }

} //namespace itensor

#endif
