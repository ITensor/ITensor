//
// Distributed under the ITensor Library License, Version 1.0.
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
    getSiP(int i) const;

    virtual IQTensor
    makeTReverse(int i) const;

    virtual IQTensor
    makeNup(int i) const;

    virtual IQTensor
    makeNdn(int i) const;

    virtual IQTensor
    makeNtot(int i) const;

    virtual IQTensor
    makeCup(int i) const;

    virtual IQTensor
    makeCdagup(int i) const;

    virtual IQTensor
    makeCdn(int i) const;

    virtual IQTensor
    makeCdagdn(int i) const;

    virtual IQTensor
    makeAup(int i) const;

    virtual IQTensor
    makeAdagup(int i) const;

    virtual IQTensor
    makeAdn(int i) const;

    virtual IQTensor
    makeAdagdn(int i) const;

    virtual IQTensor
    makeFermiPhase(int i) const;

    virtual IQTensor
    makeSz(int i) const;

    virtual IQTensor
    makeSp(int i) const;

    virtual IQTensor
    makeSm(int i) const;

    virtual IQTensor
    makeSx(int i) const;

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

inline void tJ::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        site_.at(j) = IQIndex(nameint("tJ site=",j),
            Index(nameint("Emp for site ",j),1,Site),  QN( 0,0,0),
            Index(nameint("Up for site ",j),1,Site),   QN(+1,1,1),
            Index(nameint("Dn for site ",j),1,Site),   QN(-1,1,1));
    }

inline void tJ::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

inline void tJ::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

inline int tJ::
getN() const
    { return N_; }

inline const IQIndex& tJ::
getSi(int i) const
    { return site_.at(i); }

inline IQIndex tJ::
getSiP(int i) const
    { return primed(site_.at(i)); }

inline IQIndexVal tJ::
Emp(int i) const
    {
    return getSi(i)(1);
    }

inline IQIndexVal tJ::
Up(int i) const
    {
    return getSi(i)(2);
    }

inline IQIndexVal tJ::
Dn(int i) const
    {
    return getSi(i)(3);
    }

inline IQIndexVal tJ::
EmpP(int i) const
    {
    return getSiP(i)(1);
    }

inline IQIndexVal tJ::
UpP(int i) const
    {
    return getSiP(i)(2);
    }

inline IQIndexVal tJ::
DnP(int i) const
    {
    return getSiP(i)(3);
    }

inline IQTensor tJ::
makeTReverse(int i) const
    { 
    IQTensor tr(conj(si(i)),siP(i));
    //tr(Dn(i),UpP(i)) = -1;
    tr(Dn(i),UpP(i)) = +1;
    tr(Up(i),DnP(i)) = 1;
    tr(Emp(i),EmpP(i)) = 1;
    return tr;
    }

inline IQTensor tJ::
makeNup(int i) const
    {
    IQTensor Nup(conj(si(i)),siP(i));
    Nup(Up(i),UpP(i)) = 1;
    return Nup;
    }

inline IQTensor tJ::
makeNdn(int i) const
    {
    IQTensor Ndn(conj(si(i)),siP(i));
    Ndn(Dn(i),DnP(i)) = 1;
    return Ndn;
    }

inline IQTensor tJ::
makeNtot(int i) const
    {
    IQTensor Ntot(conj(si(i)),siP(i));
    Ntot(Up(i),UpP(i)) = 1;
    Ntot(Dn(i),DnP(i)) = 1;
    return Ntot;
    }

inline IQTensor tJ::
makeCup(int i) const
    {
    IQTensor Cup(conj(si(i)),siP(i));
    Cup(Up(i),EmpP(i)) = 1;
    return Cup;
    }

inline IQTensor tJ::
makeCdagup(int i) const
    {
    IQTensor Cdagup(conj(si(i)),siP(i));
    Cdagup(Emp(i),UpP(i)) = 1;
    return Cdagup;
    }

inline IQTensor tJ::
makeCdn(int i) const
    {
    IQTensor Cdn(conj(si(i)),siP(i));
    Cdn(Dn(i),EmpP(i)) = 1;
    return Cdn;
    }

inline IQTensor tJ::
makeCdagdn(int i) const
    {
    IQTensor Cdagdn(conj(si(i)),siP(i));
    Cdagdn(Emp(i),DnP(i)) = 1;
    return Cdagdn;
    }

inline IQTensor tJ::
makeAup(int i) const
    {
    IQTensor Aup(conj(si(i)),siP(i));
    Aup(Up(i),EmpP(i)) = 1;
    return Aup;
    }

inline IQTensor tJ::
makeAdagup(int i) const
    {
    IQTensor Adagup(conj(si(i)),siP(i));
    Adagup(Emp(i),UpP(i)) = 1;
    return Adagup;
    }

inline IQTensor tJ::
makeAdn(int i) const
    {
    IQTensor Adn(conj(si(i)),siP(i));
    Adn(Dn(i),EmpP(i)) = 1;
    return Adn;
    }

inline IQTensor tJ::
makeAdagdn(int i) const
    {
    IQTensor Adagdn(conj(si(i)),siP(i));
    Adagdn(Emp(i),DnP(i)) = 1;
    return Adagdn;
    }

inline IQTensor tJ::
makeFermiPhase(int i) const
    {
    IQTensor fermiPhase(conj(si(i)),siP(i));
    fermiPhase(Emp(i),EmpP(i)) = +1;
    fermiPhase(Up(i),UpP(i)) = -1;
    fermiPhase(Dn(i),DnP(i)) = -1;
    return fermiPhase;
    }

inline IQTensor tJ::
makeSz(int i) const
    {
    IQTensor Sz(conj(si(i)),siP(i));
    Sz(Up(i),UpP(i)) = +0.5; 
    Sz(Dn(i),DnP(i)) = -0.5;
    return Sz;
    }

inline IQTensor tJ::
makeSp(int i) const
    {
    IQTensor Sp(conj(si(i)),siP(i));
    Sp(Dn(i),UpP(i)) = 1;
    return Sp;
    }

inline IQTensor tJ::
makeSm(int i) const
    {
    IQTensor Sp(conj(si(i)),siP(i));
    Sp(Up(i),DnP(i)) = 1;
    return Sp;
    }

inline IQTensor tJ::
makeSx(int i) const
    {
    IQTensor Sx(conj(si(i)),siP(i));
    Sx(Up(i),DnP(i)) = 1;
    Sx(Dn(i),UpP(i)) = 1;
    return Sx;
    }

#endif
