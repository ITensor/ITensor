//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TJ_H
#define __ITENSOR_TJ_H
#include "model.h"

class tJ : public Model
    {
    public:

    tJ();

    tJ(int N);

    tJ(std::ifstream& s) { doRead(s); }

    IQIndexVal
    Emp(int i) const;

    IQIndexVal
    Up(int i) const;

    IQIndexVal
    Dn(int i) const;

    IQIndexVal
    EmpP(int i) const;

    IQIndexVal
    UpP(int i) const;

    IQIndexVal
    DnP(int i) const;

    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndex
    getSiP(int i) const { return primed(getSi(i)); }

    IQIndexVal
    getState(int i, const String& state) const;

    IQTensor
    getOp(int i, const String& opname, const OptSet& opts) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();

    //Older interface included for backwards compatibility
    IQTensor
    makeTReverse(int i) const;
    IQTensor
    makeNup(int i) const;
    IQTensor
    makeNdn(int i) const;
    IQTensor
    makeNtot(int i) const;
    IQTensor
    makeCup(int i) const;
    IQTensor
    makeCdagup(int i) const;
    IQTensor
    makeCdn(int i) const;
    IQTensor
    makeCdagdn(int i) const;
    IQTensor
    makeAup(int i) const;
    IQTensor
    makeAdagup(int i) const;
    IQTensor
    makeAdn(int i) const;
    IQTensor
    makeAdagdn(int i) const;
    IQTensor
    makeFermiPhase(int i) const;
    IQTensor
    makeSz(int i) const;
    IQTensor
    makeSp(int i) const;
    IQTensor
    makeSm(int i) const;
    IQTensor
    makeSx(int i) const;

        
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
            Index(nameint("Emp for site ",j),1,Site),  QN( 0,0,0),
            Index(nameint("Up for site ",j),1,Site),   QN(+1,1,1),
            Index(nameint("Dn for site ",j),1,Site),   QN(-1,1,1));
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
        {
        Error("State " + state + " not recognized");
        return IQIndexVal();
        }
    }

IQIndexVal inline tJ::
Emp(int i) const
    {
    return getSi(i)(1);
    }

IQIndexVal inline tJ::
Up(int i) const
    {
    return getSi(i)(2);
    }

IQIndexVal inline tJ::
Dn(int i) const
    {
    return getSi(i)(3);
    }

IQIndexVal inline tJ::
EmpP(int i) const
    {
    return primed(getSi(i))(1);
    }

IQIndexVal inline tJ::
UpP(int i) const
    {
    return primed(getSi(i))(2);
    }

IQIndexVal inline tJ::
DnP(int i) const
    {
    return primed(getSiP(i))(3);
    }

IQTensor inline tJ::
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
               DnP(sP(3));

    IQTensor Op(conj(s),sP);

    if(opname == "TReverse")
        {
        Op(Em,EmP) = +1;
        Op(Dn,UpP) = -1; //correct?
        Op(Up,DnP) = +1;
        }
    else
    if(opname == "Nup")
        {
        Op(Up,UpP) = 1;
        }
    else
    if(opname == "Ndn")
        {
        Op(Dn,DnP) = 1;
        }
    else
    if(opname == "Ntot")
        {
        Op(Up,UpP) = 1;
        Op(Dn,DnP) = 1;
        }
    else
    if(opname == "Cup")
        {
        Op(Up,EmP) = 1; 
        }
    else
    if(opname == "Cdagup")
        {
        Op(Em,UpP) = 1; 
        }
    else
    if(opname == "Cdn")
        {
        Op(Dn,EmP) = 1; 
        }
    else
    if(opname == "Cdagdn")
        {
        Op(Em,DnP) = 1; 
        }
    else
    if(opname == "Aup")
        {
        Op(Up,EmP) = 1; 
        }
    else
    if(opname == "Adagup")
        {
        Op(Em,UpP) = 1; 
        }
    else
    if(opname == "Adn")
        {
        Op(Dn,EmP) = 1; 
        }
    else
    if(opname == "Adagdn")
        {
        Op(Em,DnP) = 1; 
        }
    else
    if(opname == "FermiPhase" || opname == "F")
        {
        Op(Em,EmP) = +1; 
        Op(Up,UpP) = -1;
        Op(Dn,DnP) = -1;
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
    if(opname == "Sp")
        {
        Op(Dn,UpP) = 1;
        }
    else
    if(opname == "Sm")
        {
        Op(Up,DnP) = 1;
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }


//
// Older interface (deprecated)
// included for backwards compatibility
//

IQTensor inline tJ::
makeTReverse(int i) const
    { 
    IQTensor tr(conj(si(i)),siP(i));
    //tr(Dn(i),UpP(i)) = -1;
    tr(Dn(i),UpP(i)) = +1;
    tr(Up(i),DnP(i)) = 1;
    tr(Emp(i),EmpP(i)) = 1;
    return tr;
    }

IQTensor inline tJ::
makeNup(int i) const
    {
    IQTensor Nup(conj(si(i)),siP(i));
    Nup(Up(i),UpP(i)) = 1;
    return Nup;
    }

IQTensor inline tJ::
makeNdn(int i) const
    {
    IQTensor Ndn(conj(si(i)),siP(i));
    Ndn(Dn(i),DnP(i)) = 1;
    return Ndn;
    }

IQTensor inline tJ::
makeNtot(int i) const
    {
    IQTensor Ntot(conj(si(i)),siP(i));
    Ntot(Up(i),UpP(i)) = 1;
    Ntot(Dn(i),DnP(i)) = 1;
    return Ntot;
    }

IQTensor inline tJ::
makeCup(int i) const
    {
    IQTensor Cup(conj(si(i)),siP(i));
    Cup(Up(i),EmpP(i)) = 1;
    return Cup;
    }

IQTensor inline tJ::
makeCdagup(int i) const
    {
    IQTensor Cdagup(conj(si(i)),siP(i));
    Cdagup(Emp(i),UpP(i)) = 1;
    return Cdagup;
    }

IQTensor inline tJ::
makeCdn(int i) const
    {
    IQTensor Cdn(conj(si(i)),siP(i));
    Cdn(Dn(i),EmpP(i)) = 1;
    return Cdn;
    }

IQTensor inline tJ::
makeCdagdn(int i) const
    {
    IQTensor Cdagdn(conj(si(i)),siP(i));
    Cdagdn(Emp(i),DnP(i)) = 1;
    return Cdagdn;
    }

IQTensor inline tJ::
makeAup(int i) const
    {
    IQTensor Aup(conj(si(i)),siP(i));
    Aup(Up(i),EmpP(i)) = 1;
    return Aup;
    }

IQTensor inline tJ::
makeAdagup(int i) const
    {
    IQTensor Adagup(conj(si(i)),siP(i));
    Adagup(Emp(i),UpP(i)) = 1;
    return Adagup;
    }

IQTensor inline tJ::
makeAdn(int i) const
    {
    IQTensor Adn(conj(si(i)),siP(i));
    Adn(Dn(i),EmpP(i)) = 1;
    return Adn;
    }

IQTensor inline tJ::
makeAdagdn(int i) const
    {
    IQTensor Adagdn(conj(si(i)),siP(i));
    Adagdn(Emp(i),DnP(i)) = 1;
    return Adagdn;
    }

IQTensor inline tJ::
makeFermiPhase(int i) const
    {
    IQTensor fermiPhase(conj(si(i)),siP(i));
    fermiPhase(Emp(i),EmpP(i)) = +1;
    fermiPhase(Up(i),UpP(i)) = -1;
    fermiPhase(Dn(i),DnP(i)) = -1;
    return fermiPhase;
    }

IQTensor inline tJ::
makeSz(int i) const
    {
    IQTensor Sz(conj(si(i)),siP(i));
    Sz(Up(i),UpP(i)) = +0.5; 
    Sz(Dn(i),DnP(i)) = -0.5;
    return Sz;
    }

IQTensor inline tJ::
makeSp(int i) const
    {
    IQTensor Sp(conj(si(i)),siP(i));
    Sp(Dn(i),UpP(i)) = 1;
    return Sp;
    }

IQTensor inline tJ::
makeSm(int i) const
    {
    IQTensor Sp(conj(si(i)),siP(i));
    Sp(Up(i),DnP(i)) = 1;
    return Sp;
    }

IQTensor inline tJ::
makeSx(int i) const
    {
    IQTensor Sx(conj(si(i)),siP(i));
    Sx(Up(i),DnP(i)) = 1;
    Sx(Dn(i),UpP(i)) = 1;
    return Sx;
    }

#endif
