//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINONE_H
#define __ITENSOR_SPINONE_H
#include "../model.h"

class SpinOne : public Model
    {
    public:

    SpinOne();

    SpinOne(int N, bool shalf_edge = false);

    SpinOne(std::ifstream& s) { doRead(s); }

    IQIndexVal
    Up(int i) const;

    IQIndexVal
    Z0(int i) const;

    IQIndexVal
    Dn(int i) const;

    IQIndexVal
    UpP(int i) const;

    IQIndexVal
    Z0P(int i) const;

    IQIndexVal
    DnP(int i) const;

    //(Sz)^2, etc. operators

    IQTensor
    sz2(int i) const { return makeSx2(i); }

    IQTensor
    sx2(int i) const { return makeSx2(i); }

    IQTensor
    sy2(int i) const { return makeSx2(i); }

    private:

    virtual int
    getNN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndex
    getSiP(int i) const;

    virtual IQTensor
    makeSz(int i) const;

    virtual IQTensor
    makeSx(int i) const;

    virtual IQTensor
    makeISy(int i) const;

    virtual IQTensor
    makeSp(int i) const;

    virtual IQTensor
    makeSm(int i) const;

    virtual IQTensor
    makeSz2(int i) const;

    virtual IQTensor
    makeSx2(int i) const;

    virtual IQTensor
    makeSy2(int i) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites(bool shalf_edge);
        
    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;

    };

inline SpinOne::
SpinOne()
    : N_(-1)
    { }

inline SpinOne::
SpinOne(int N, bool shalf_edge)
    : N_(N),
      site_(N_+1)
    { 
    constructSites(shalf_edge);
    }

inline void SpinOne::
constructSites(bool shalf_edge)
    {
    for(int j = 1; j <= N_; ++j)
        if(shalf_edge && (j == 1 || j == N_))
            {
            site_.at(j) = IQIndex(nameint("S=1/2, site=",j),
                Index(nameint("Up for site",j),1,Site),QN(+1,0),
                Index(nameint("Dn for site",j),1,Site),QN(-1,0));
            }
        else
            {
            site_.at(j) = IQIndex(nameint("S=1, site=",j),
                Index(nameint("Up for site",j),1,Site),QN(+2,0),
                Index(nameint("Z0 for site",j),1,Site),QN( 0,0),
                Index(nameint("Dn for site",j),1,Site),QN(-2,0));
            }
    }

inline void SpinOne::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

inline void SpinOne::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

inline int SpinOne::
getNN() const
    { return N_; }

inline const IQIndex& SpinOne::
getSi(int i) const
    { return site_.at(i); }

inline IQIndex SpinOne::
getSiP(int i) const
    { return site_.at(i).primed(); }

inline IQIndexVal SpinOne::
Up(int i) const
    {
    return getSi(i)(1);
    }

inline IQIndexVal SpinOne::
Z0(int i) const
    {
    IQIndex ind = getSi(i);
    if(ind.m() == 2)
        Error("Z0 not defined for spin 1/2 site");
    return ind(2);
    }

inline IQIndexVal SpinOne::
Dn(int i) const
    {
    IQIndex ind = getSi(i);
    return ind(ind.m());
    }

inline IQIndexVal SpinOne::
UpP(int i) const
    {
    return getSiP(i)(1);
    }

inline IQIndexVal SpinOne::
Z0P(int i) const
    {
    IQIndex ind = getSiP(i);
    if(ind.m() == 2)
        Error("Z0 not defined for spin 1/2 site");
    return ind(2);
    }

inline IQIndexVal SpinOne::
DnP(int i) const
    {
    IQIndex ind = getSiP(i);
    return ind(ind.m());
    }

inline IQTensor SpinOne::
makeSz(int i) const
    {
    IQTensor Sz(conj(si(i)),siP(i));
    if(si(i).m() == 2)
        {
        Sz(Up(i),UpP(i)) = +0.5;
        Sz(Dn(i),DnP(i)) = -0.5;
        }
    else
        {
        Sz(Up(i),UpP(i)) = +1.;
        Sz(Dn(i),DnP(i)) = -1.;
        }
    return Sz;
    }

inline IQTensor SpinOne::
makeSx(int i) const
    {
    IQTensor Sx(conj(si(i)),siP(i));
    if(si(i).m() == 2)
        {
        Sx(Up(i),DnP(i)) = +0.5;
        Sx(Dn(i),UpP(i)) = +0.5;
        }
    else
        {
        Sx(Up(i),Z0P(i)) = ISqrt2; 
        Sx(Z0(i),UpP(i)) = ISqrt2;
        Sx(Z0(i),DnP(i)) = ISqrt2; 
        Sx(Dn(i),Z0P(i)) = ISqrt2;
        }
    return Sx;
    }

inline IQTensor SpinOne::
makeISy(int i) const
    {
    IQTensor ISy(conj(si(i)),siP(i));
    if(si(i).m() == 2)
        {
        ISy(Up(i),DnP(i)) = -0.5;
        ISy(Dn(i),UpP(i)) = +0.5;
        }
    else
        {
        ISy(Up(i),Z0P(i)) = +ISqrt2; 
        ISy(Z0(i),UpP(i)) = -ISqrt2;
        ISy(Z0(i),DnP(i)) = +ISqrt2; 
        ISy(Dn(i),Z0P(i)) = -ISqrt2;
        }
    return ISy;
    }

inline IQTensor SpinOne::
makeSp(int i) const
    {
    IQTensor Sp(conj(si(i)),siP(i));
    if(si(i).m() == 2)
        {
        Sp(Dn(i),UpP(i)) = 1;
        }
    else
        {
        Sp(Dn(i),Z0P(i)) = Sqrt2; 
        Sp(Z0(i),UpP(i)) = Sqrt2;
        }
    return Sp;
    }

inline IQTensor SpinOne::
makeSm(int i) const
    {
    IQTensor Sm(conj(si(i)),siP(i));
    if(si(i).m() == 2)
        {
        Sm(Up(i),DnP(i)) = 1;
        }
    else
        {
        Sm(Up(i),Z0P(i)) = Sqrt2; 
        Sm(Z0(i),DnP(i)) = Sqrt2;
        }
    return Sm;
    }


inline IQTensor SpinOne::
makeSz2(int i) const
    {
    if(si(i).m() == 2) Error("Sz^2 only non-trivial for S=1 sites");
    IQTensor Sz2(conj(si(i)),siP(i));
    Sz2(Up(i),UpP(i)) = 1; 
    Sz2(Dn(i),DnP(i)) = 1;
    return Sz2;
    }

inline IQTensor SpinOne::
makeSx2(int i) const
    {
    if(si(i).m() == 2) Error("Sx^2 only non-trivial for S=1 sites");
    IQTensor Sx2(conj(si(i)),siP(i));
    Sx2(Up(i),UpP(i)) = 0.5; 
    Sx2(Up(i),DnP(i)) = 0.5;
    Sx2(Z0(i),Z0P(i)) = 1;
    Sx2(Dn(i),DnP(i)) = 0.5; 
    Sx2(Dn(i),UpP(i)) = 0.5;
    return Sx2;
    }

inline IQTensor SpinOne::
makeSy2(int i) const
    {
    if(si(i).m() == 2) Error("Sy^2 only non-trivial for S=1 sites");
    IQTensor Sy2(conj(si(i)),siP(i));
    Sy2(Up(i),UpP(i)) = 0.5; 
    Sy2(Up(i),DnP(i)) = -0.5;
    Sy2(Z0(i),Z0P(i)) = 1;
    Sy2(Dn(i),DnP(i)) = 0.5; 
    Sy2(Dn(i),UpP(i)) = -0.5;
    return Sy2;
    }

#endif
