//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HUBBARD_H
#define __ITENSOR_HUBBARD_H
#include "../model.h"

class Hubbard : public Model
    {
    public:

    Hubbard();

    Hubbard(int N, 
            const OptSet& opts = Global::opts());

    bool
    conserveNf() const { return conserveNf_; }

    IQIndexVal
    Emp(int i) const;

    IQIndexVal
    Up(int i) const;

    IQIndexVal
    Dn(int i) const;

    IQIndexVal
    UpDn(int i) const;

    IQIndexVal
    EmpP(int i) const;

    IQIndexVal
    UpP(int i) const;

    IQIndexVal
    DnP(int i) const;

    IQIndexVal
    UpDnP(int i) const;

    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const OptSet& opts = Global::opts()) const;

    IQTensor
    makeTReverse(int i) const { return getOp(i,"TReverse"); }

    IQTensor
    makeNup(int i) const { return getOp(i,"Nup"); }

    IQTensor
    makeNdn(int i) const { return getOp(i,"Ndn"); }

    IQTensor
    makeNupdn(int i) const { return getOp(i,"Nupdn"); }

    IQTensor
    makeNtot(int i) const { return getOp(i,"Ntot"); }

    IQTensor
    makeCup(int i) const { return getOp(i,"Cup"); }

    IQTensor
    makeCdagup(int i) const { return getOp(i,"Cdagup"); }

    IQTensor
    makeCdn(int i) const { return getOp(i,"Cdn"); }

    IQTensor
    makeCdagdn(int i) const { return getOp(i,"Cdagdn"); }

    IQTensor
    makeAup(int i) const { return getOp(i,"Aup"); }

    IQTensor
    makeAdagup(int i) const { return getOp(i,"Adagup"); }

    IQTensor
    makeAdn(int i) const { return getOp(i,"Adn"); }

    IQTensor
    makeAdagdn(int i) const { return getOp(i,"Adagdn"); }

    IQTensor
    makeFermiPhase(int i) const { return getOp(i,"F"); }

    IQTensor
    makeSz(int i) const { return getOp(i,"Sz"); }

    IQTensor
    makeSx(int i) const { return getOp(i,"Sx"); }

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();
        
    //Data members -----------------

    int N_;
    bool conserveNf_;

    std::vector<IQIndex> site_;

    };

inline Hubbard::
Hubbard()
    : N_(-1),
    conserveNf_(true)
    { }

inline Hubbard::
Hubbard(int N, const OptSet& opts)
    : N_(N),
      site_(N_+1)
    { 
    conserveNf_ = opts.getBool("ConserveNf",true);
    constructSites();
    }

void inline Hubbard::
constructSites()
    {
    const int One = (conserveNf_ ? 1 : 0);
    const int Two = 2*One;
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("Hubbard site=",j),
            Index(nameint("Emp for site ",j),1,Site),  QN( 0,0,0),
            Index(nameint("Up for site ",j),1,Site),   QN(+1,One,1),
            Index(nameint("Dn for site ",j),1,Site),   QN(-1,One,1),
            Index(nameint("UpDn for site ",j),1,Site), QN( 0,Two,0));
        }
    }

void inline Hubbard::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    if(site_.at(1).qn(2).Nf() == 1)
        conserveNf_ = true;
    else
        conserveNf_ = false;
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
getState(int i, const String& state) const
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
    }

IQIndexVal inline Hubbard::
Emp(int i) const
    {
    return getState(i,"0");
    }

IQIndexVal inline Hubbard::
Up(int i) const
    {
    return getState(i,"+");
    }

IQIndexVal inline Hubbard::
Dn(int i) const
    {
    return getState(i,"-");
    }

IQIndexVal inline Hubbard::
UpDn(int i) const
    {
    return getState(i,"S");
    }

IQIndexVal inline Hubbard::
EmpP(int i) const
    {
    return primed(getState(i,"0"));
    }

IQIndexVal inline Hubbard::
UpP(int i) const
    {
    return primed(getState(i,"+"));
    }

IQIndexVal inline Hubbard::
DnP(int i) const
    {
    return primed(getState(i,"-"));
    }

IQIndexVal inline Hubbard::
UpDnP(int i) const
    {
    return primed(getState(i,"S"));
    }

inline IQTensor Hubbard::
getOp(int i, const String& opname, const OptSet& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = primed(s);

    IQIndexVal Em(s(1)),
               EmP(sP(1)),
               Up(s(2)),
               UpP(sP(2)),
               Dn(s(3)),
               DnP(sP(3)),
               UD(s(4)),
               UDP(sP(4));

    IQTensor Op(conj(s),sP);

    if(opname == "TReverse")
        {
        Op(Em,EmP) = +1;
        Op(UD,UDP) = -1;
        //should one of the following two lines be a -1?
        Op(Dn,UpP) = +1;
        Op(Up,DnP) = +1;
        }
    else
    if(opname == "Nup")
        {
        Op(Up,UpP) = 1;
        Op(UD,UDP) = 1;
        }
    else
    if(opname == "Ndn")
        {
        Op(Dn,DnP) = 1;
        Op(UD,UDP) = 1;
        }
    else
    if(opname == "Nupdn")
        {
        Op(UD,UDP) = 1;
        }
    else
    if(opname == "Ntot")
        {
        Op(Up,UpP) = 1;
        Op(Dn,DnP) = 1;
        Op(UD,UDP) = 2;
        }
    else
    if(opname == "Cup")
        {
        Op(Up,EmP) = 1; 
        Op(UD,DnP) = 1; 
        }
    else
    if(opname == "Cdagup")
        {
        Op(Em,UpP) = 1; 
        Op(Dn,UDP) = 1;
        }
    else
    if(opname == "Cdn")
        {
        Op(Dn,EmP) = 1; 
        Op(UD,UpP) = -1; 
        }
    else
    if(opname == "Cdagdn")
        {
        Op(Em,DnP) = 1; 
        Op(Up,UDP) = -1;
        }
    else
    if(opname == "Aup")
        {
        Op(Up,EmP) = 1; 
        Op(UD,DnP) = 1; 
        }
    else
    if(opname == "Adagup")
        {
        Op(Em,UpP) = 1; 
        Op(Dn,UDP) = 1;
        }
    else
    if(opname == "Adn")
        {
        Op(Dn,EmP) = 1; 
        Op(UD,UpP) = 1; 
        }
    else
    if(opname == "Adagdn")
        {
        Op(Em,DnP) = 1; 
        Op(Up,UDP) = 1;
        }
    else
    if(opname == "FermiPhase" || opname == "F")
        {
        Op(Em,EmP) = +1; 
        Op(Up,UpP) = -1;
        Op(Dn,DnP) = -1;
        Op(UD,UDP) = +1;
        }
    else
    if(opname == "Sz")
        {
        Op(Up,UpP) = +0.5; 
        Op(Dn,DnP) = -0.5;
        }
    else
    if(opname == "Sx")
        {
        Op(Up,DnP) = 1; 
        Op(Dn,UpP) = 1;
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

#endif
