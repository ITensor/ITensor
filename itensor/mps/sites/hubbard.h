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

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, String const& state) const;

    virtual IQTensor
    getOp(int i, String const& opname, Args const& args = Global::args()) const;

    DefaultOpsT
    getDefaultOps(Args const& args) const;

    void 
    constructSites();

    void
    doRead(std::istream& s);

    void
    doWrite(std::ostream& s) const;

        
    //Data members -----------------

    int N_;
    bool conserveNf_,
         conserveSz_;

    std::vector<IQIndex> site_;

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
    : N_(-1),
    conserveNf_(true),
    conserveSz_(true)
    { }

inline Hubbard::
Hubbard(int N, Args const& args)
    : N_(N),
      site_(N_+1)
    { 
    conserveNf_ = args.getBool("ConserveNf",true);
    conserveSz_ = args.getBool("ConserveSz",true);
    constructSites();
    }

void inline Hubbard::
constructSites()
    {
    int Up = (conserveSz_ ? +1 : 0),
        Dn = -Up;
    if(conserveNf_)
        {
        for(auto j : range1(N_))
            {
            site_.at(j) = IQIndex(nameint("site=",j),
                Index(nameint("Emp ",j),1,Site), QN("Sz=", 0,"Nf=",0),
                Index(nameint("Up ",j),1,Site),  QN("Sz=",Up,"Nf=",1),
                Index(nameint("Dn ",j),1,Site),  QN("Sz=",Dn,"Nf=",1),
                Index(nameint("UpDn ",j),1,Site),QN("Sz=", 0,"Nf=",2));
            }
        }
    else //don't conserve Nf, only fermion parity
        {
        if(!conserveSz_) Error("One of ConserveSz or ConserveNf must be true for Hubbard sites");

        for(auto j : range1(N_))
            {
            site_.at(j) = IQIndex(nameint("site=",j),
                Index(nameint("Emp ",j),1,Site), QN("Sz=", 0,"Pf=",0),
                Index(nameint("Up ",j),1,Site),  QN("Sz=",+1,"Pf=",1),
                Index(nameint("Dn ",j),1,Site),  QN("Sz=",-1,"Pf=",1),
                Index(nameint("UpDn ",j),1,Site),QN("Sz=", 0,"Pf=",0));
            }
        }
    }

void inline Hubbard::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);

    conserveNf_ = (site_.at(1).qn(2)(2) == 1);
    }

void inline Hubbard::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

int inline Hubbard::
getN() const
    { return N_; }

inline const IQIndex& Hubbard::
getSi(int i) const
    { return site_.at(i); }

inline IQIndexVal Hubbard::
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
    if(state == "S" || state == "UpDn") 
        {
        return getSi(i)(4);
        }
    else
        {
        Error("State " + state + " not recognized");
        return getSi(i)(1);
        }
    return IQIndexVal();
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
    static const std::vector<String> dops_(initDefaultOps());
    return dops_;
    }

} //namespace itensor

#endif
