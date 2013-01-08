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

    Hubbard(std::ifstream& s) { doRead(s); }

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

    virtual IQIndex
    getSiP(int i) const;

    virtual IQTensor
    makeTReverse(int i) const;

    virtual IQTensor
    makeNup(int i) const;

    virtual IQTensor
    makeNdn(int i) const;

    virtual IQTensor
    makeNupdn(int i) const;

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
    makeSx(int i) const;

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

IQIndex inline Hubbard::
getSiP(int i) const
    { return site_.at(i).primed(); }

IQIndexVal inline Hubbard::
Emp(int i) const
    {
    return getSi(i)(1);
    }

IQIndexVal inline Hubbard::
Up(int i) const
    {
    return getSi(i)(2);
    }

IQIndexVal inline Hubbard::
Dn(int i) const
    {
    return getSi(i)(3);
    }

IQIndexVal inline Hubbard::
UpDn(int i) const
    {
    return getSi(i)(4);
    }

IQIndexVal inline Hubbard::
EmpP(int i) const
    {
    return getSiP(i)(1);
    }

IQIndexVal inline Hubbard::
UpP(int i) const
    {
    return getSiP(i)(2);
    }

IQIndexVal inline Hubbard::
DnP(int i) const
    {
    return getSiP(i)(3);
    }

IQIndexVal inline Hubbard::
UpDnP(int i) const
    {
    return getSiP(i)(4);
    }

IQTensor inline Hubbard::
makeTReverse(int i) const
    { 
    IQTensor tr(conj(si(i)),siP(i));
    tr(UpDn(i),UpDnP(i)) = -1;
    //tr(Dn(i),UpP(i)) = -1;
    tr(Dn(i),UpP(i)) = +1;
    tr(Up(i),DnP(i)) = 1;
    tr(Emp(i),EmpP(i)) = 1;
    return tr;
    }

IQTensor inline Hubbard::
makeNup(int i) const
    {
    IQTensor Nup(conj(si(i)),siP(i));
    Nup(Up(i),UpP(i)) = 1;
    Nup(UpDn(i),UpDnP(i)) = 1;
    return Nup;
    }

IQTensor inline Hubbard::
makeNdn(int i) const
    {
    IQTensor Ndn(conj(si(i)),siP(i));
    Ndn(Dn(i),DnP(i)) = 1;
    Ndn(UpDn(i),UpDnP(i)) = 1;
    return Ndn;
    }

IQTensor inline Hubbard::
makeNupdn(int i) const
    {
    IQTensor Nupdn(conj(si(i)),siP(i));
    Nupdn(UpDn(i),UpDnP(i)) = 1;
    return Nupdn;
    }

IQTensor inline Hubbard::
makeNtot(int i) const
    {
    IQTensor Ntot(conj(si(i)),siP(i));
    Ntot(Up(i),UpP(i)) = 1;
    Ntot(Dn(i),DnP(i)) = 1;
    Ntot(UpDn(i),UpDnP(i)) = 2;
    return Ntot;
    }

IQTensor inline Hubbard::
makeCup(int i) const
    {
    IQTensor Cup(conj(si(i)),siP(i));
    Cup(Up(i),EmpP(i)) = 1;
    Cup(UpDn(i),DnP(i)) = -1;
    return Cup;
    }

IQTensor inline Hubbard::
makeCdagup(int i) const
    {
    IQTensor Cdagup(conj(si(i)),siP(i));
    Cdagup(Emp(i),UpP(i)) = 1;
    Cdagup(Dn(i),UpDnP(i)) = -1;
    return Cdagup;
    }

IQTensor inline Hubbard::
makeCdn(int i) const
    {
    IQTensor Cdn(conj(si(i)),siP(i));
    Cdn(Dn(i),EmpP(i)) = 1;
    Cdn(UpDn(i),UpP(i)) = 1;
    return Cdn;
    }

IQTensor inline Hubbard::
makeCdagdn(int i) const
    {
    IQTensor Cdagdn(conj(si(i)),siP(i));
    Cdagdn(Emp(i),DnP(i)) = 1;
    Cdagdn(Up(i),UpDnP(i)) = 1;
    return Cdagdn;
    }

IQTensor inline Hubbard::
makeAup(int i) const
    {
    IQTensor Aup(conj(si(i)),siP(i));
    Aup(Up(i),EmpP(i)) = 1;
    Aup(UpDn(i),DnP(i)) = 1;
    return Aup;
    }

IQTensor inline Hubbard::
makeAdagup(int i) const
    {
    IQTensor Adagup(conj(si(i)),siP(i));
    Adagup(Emp(i),UpP(i)) = 1;
    Adagup(Dn(i),UpDnP(i)) = 1;
    return Adagup;
    }

IQTensor inline Hubbard::
makeAdn(int i) const
    {
    IQTensor Adn(conj(si(i)),siP(i));
    Adn(Dn(i),EmpP(i)) = 1;
    Adn(UpDn(i),UpP(i)) = 1;
    return Adn;
    }

IQTensor inline Hubbard::
makeAdagdn(int i) const
    {
    IQTensor Adagdn(conj(si(i)),siP(i));
    Adagdn(Emp(i),DnP(i)) = 1;
    Adagdn(Up(i),UpDnP(i)) = 1;
    return Adagdn;
    }

IQTensor inline Hubbard::
makeFermiPhase(int i) const
    {
    IQTensor fermiPhase(conj(si(i)),siP(i));
    fermiPhase(Emp(i),EmpP(i)) = +1;
    fermiPhase(Up(i),UpP(i)) = -1;
    fermiPhase(Dn(i),DnP(i)) = -1;
    fermiPhase(UpDn(i),UpDnP(i)) = +1;
    return fermiPhase;
    }

IQTensor inline Hubbard::
makeSz(int i) const
    {
    IQTensor Sz(conj(si(i)),siP(i));
    Sz(Up(i),UpP(i)) = +0.5; 
    Sz(Dn(i),DnP(i)) = -0.5;
    return Sz;
    }

IQTensor inline Hubbard::
makeSx(int i) const
    {
    IQTensor Sx(conj(si(i)),siP(i));
    Sx(Up(i),DnP(i)) = 1;
    Sx(Dn(i),UpP(i)) = 1;
    return Sx;
    }

#endif
