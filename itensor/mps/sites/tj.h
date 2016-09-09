//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TJ_H
#define __ITENSOR_TJ_H
#include "itensor/mps/siteset.h"

namespace itensor {

class tJ : public SiteSet
    {
    public:

    tJ();

    tJ(int N);

    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndex
    getSiP(int i) const { return prime(getSi(i)); }

    IQIndexVal
    getState(int i, String const& state) const;

    IQTensor
    getOp(int i, String const& opname, Args const& args) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();

    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;

    };

inline tJ::
tJ()
    : N_(-1)
    { }

inline tJ::
tJ(int N)
    : N_(N),
      site_(N_+1)
    { 
    constructSites();
    }

void inline tJ::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        site_.at(j) = IQIndex(nameint("tJ site=",j),
            Index(nameint("Emp ",j),1,Site), QN("Sz=", 0,"Nf=",0),
            Index(nameint("Up ",j),1,Site),  QN("Sz=",+1,"Nf=",1),
            Index(nameint("Dn ",j),1,Site),  QN("Sz=",-1,"Nf=",1));
    }

void inline tJ::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

void inline tJ::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

int inline tJ::
getN() const
    { return N_; }

inline const IQIndex& tJ::
getSi(int i) const
    { return site_.at(i); }

IQIndexVal inline tJ::
getState(int i, String const& state) const
    {
    if(state == "0" || state == "Emp") 
        {
        return getSi(i)(1);
        }
    else 
    if(state == "+" || state == "Up") 
        {
        return getSi(i)(2);
        }
    else 
    if(state == "-" || state == "Dn") 
        {
        return getSi(i)(3);
        }
    else
        {
        Error("State " + state + " not recognized");
        return IQIndexVal();
        }
    }


IQTensor inline tJ::
getOp(int i, String const& opname, Args const& args) const
    {
    auto s = si(i);
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
    if(opname == "S+" || opname == "Sp")
        {
        Op.set(Dn,UpP,1);
        }
    else
    if(opname == "S-" || opname == "Sm")
        {
        Op.set(Up,DnP,1);
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

} //namespace itensor

#endif
