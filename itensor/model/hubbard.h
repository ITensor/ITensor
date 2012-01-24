#ifndef __ITENSOR_HUBBARD_H
#define __ITENSOR_HUBBARD_H
#include "../model.h"

class Hubbard : public Model
    {
    public:

    Hubbard();

    Hubbard(int N);

    Hubbard(std::ifstream& s) { doRead(s); }

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
    getNN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndex
    getSiP(int i) const;

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
    makeFermiPhase(int i) const;

    virtual IQTensor
    makeSz(int i) const;

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

inline Hubbard::
Hubbard()
    : N_(-1)
    { }

inline Hubbard::
Hubbard(int N)
    : N_(N),
      site_(N_+1)
    { 
    constructSites();
    }

inline void Hubbard::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("Hubbard, site=",i),
            Index(nameint("Emp for site ",i),1,Site),  QN( 0,0,0),
            Index(nameint("Up for site ",i),1,Site),   QN(+1,1,1),
            Index(nameint("Dn for site ",i),1,Site),   QN(-1,1,1),
            Index(nameint("UpDn for site ",i),1,Site), QN( 0,2,0));
        }
    }

inline void Hubbard::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

inline void Hubbard::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

inline int Hubbard::
getNN() const
    { return N_; }

inline const IQIndex& Hubbard::
getSi(int i) const
    { return site_.at(i); }

inline IQIndex Hubbard::
getSiP(int i) const
    { return site_.at(i).primed(); }

inline IQIndexVal Hubbard::
Emp(int i) const
    {
    return getSi(i)(1);
    }

inline IQIndexVal Hubbard::
Up(int i) const
    {
    return getSi(i)(2);
    }

inline IQIndexVal Hubbard::
Dn(int i) const
    {
    return getSi(i)(3);
    }

inline IQIndexVal Hubbard::
UpDn(int i) const
    {
    return getSi(i)(4);
    }

inline IQIndexVal Hubbard::
EmpP(int i) const
    {
    return getSiP(i)(1);
    }

inline IQIndexVal Hubbard::
UpP(int i) const
    {
    return getSiP(i)(2);
    }

inline IQIndexVal Hubbard::
DnP(int i) const
    {
    return getSiP(i)(3);
    }

inline IQIndexVal Hubbard::
UpDnP(int i) const
    {
    return getSiP(i)(4);
    }

inline IQTensor Hubbard::
makeNup(int i) const
    {
    IQTensor Nup(conj(si(i)),siP(i));
    Nup(Up(i),UpP(i)) = 1;
    Nup(UpDn(i),UpDnP(i)) = 1;
    return Nup;
    }

inline IQTensor Hubbard::
makeNdn(int i) const
    {
    IQTensor Ndn(conj(si(i)),siP(i));
    Ndn(Dn(i),DnP(i)) = 1;
    Ndn(UpDn(i),UpDnP(i)) = 1;
    return Ndn;
    }

inline IQTensor Hubbard::
makeNupdn(int i) const
    {
    IQTensor Nupdn(conj(si(i)),siP(i));
    Nupdn(UpDn(i),UpDnP(i)) = 1;
    return Nupdn;
    }

inline IQTensor Hubbard::
makeNtot(int i) const
    {
    IQTensor Ntot(conj(si(i)),siP(i));
    Ntot(Up(i),UpP(i)) = 1;
    Ntot(Dn(i),DnP(i)) = 1;
    Ntot(UpDn(i),UpDnP(i)) = 2;
    return Ntot;
    }

inline IQTensor Hubbard::
makeCup(int i) const
    {
    IQTensor Cup(conj(si(i)),siP(i));
    Cup(Up(i),EmpP(i)) = 1;
    Cup(UpDn(i),DnP(i)) = -1;
    return Cup;
    }

inline IQTensor Hubbard::
makeCdagup(int i) const
    {
    IQTensor Cdagup(conj(si(i)),siP(i));
    Cdagup(Emp(i),UpP(i)) = 1;
    Cdagup(Dn(i),UpDnP(i)) = -1;
    return Cdagup;
    }

inline IQTensor Hubbard::
makeCdn(int i) const
    {
    IQTensor Cdn(conj(si(i)),siP(i));
    Cdn(Dn(i),EmpP(i)) = 1;
    Cdn(UpDn(i),UpP(i)) = 1;
    return Cdn;
    }

inline IQTensor Hubbard::
makeCdagdn(int i) const
    {
    IQTensor Cdagdn(conj(si(i)),siP(i));
    Cdagdn(Emp(i),DnP(i)) = 1;
    Cdagdn(Up(i),UpDnP(i)) = 1;
    return Cdagdn;
    }

inline IQTensor Hubbard::
makeFermiPhase(int i) const
    {
    IQTensor fermiPhase(conj(si(i)),siP(i));
    fermiPhase(Emp(i),EmpP(i)) = +1;
    fermiPhase(Up(i),UpP(i)) = -1;
    fermiPhase(Dn(i),DnP(i)) = -1;
    fermiPhase(UpDn(i),UpDnP(i)) = +1;
    return fermiPhase;
    }

inline IQTensor Hubbard::
makeSz(int i) const
    {
    IQTensor Sz(conj(si(i)),siP(i));
    Sz(Up(i),UpP(i)) = +0.5; 
    Sz(Dn(i),DnP(i)) = -0.5;
    return Sz;
    }

#endif
