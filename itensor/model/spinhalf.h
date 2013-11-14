//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINHALF_H
#define __ITENSOR_SPINHALF_H
#include "../model.h"

class SpinHalf : public Model
    {
    public:

    SpinHalf();

    SpinHalf(int N);

    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const OptSet& opts) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();

    //Older interface included for backwards compatibility
    public:
    IQIndexVal
    Up(int i) const;
    IQIndexVal
    Dn(int i) const;
    IQIndexVal
    UpP(int i) const;
    IQIndexVal
    DnP(int i) const;
    private:
    IQTensor
    makeSz(int i) const;
    IQTensor
    makeSx(int i) const;
    IQTensor
    makeISy(int i) const;
    IQTensor
    makeSp(int i) const;
    IQTensor
    makeSm(int i) const;

    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;
        
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
        site_.at(j) = IQIndex(nameint("S=1/2 site=",j),
            Index(nameint("Up for site",j),1,Site),QN(+1,0),
            Index(nameint("Dn for site",j),1,Site),QN(-1,0));
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
getOp(int i, const String& opname, const OptSet& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = primed(s);

    IQIndexVal Up(s(1)),
               UpP(sP(1)),
               Dn(s(2)),
               DnP(sP(2));

    IQTensor Op(conj(s),sP);

    if(opname == "Sz")
        {
        Op(Up,UpP) = +0.5;
        Op(Dn,DnP) = -0.5;
        }
    else
    if(opname == "Sx")
        {
        Op(Up,DnP) = +0.5;
        Op(Dn,UpP) = +0.5;
        }
    else
    if(opname == "ISy")
        {
        Op(Up,DnP) = -0.5;
        Op(Dn,UpP) = +0.5;
        }
    else
    if(opname == "Sy")
        {
        Op(Up,DnP) = +0.5;
        Op(Dn,UpP) = -0.5;
        Op *= Complex_i;
        }
    else
    if(opname == "Sp" || opname == "S-")
        {
        Op(Dn,UpP) = 1;
        }
    else
    if(opname == "Sm" || opname == "S+")
        {
        Op(Up,DnP) = 1;
        }
    else
    if(opname == "projUp")
        {
        Op(Up,UpP) = 1; 
        }
    else
    if(opname == "projDn")
        {
        Op(Dn,DnP) = 1; 
        }
    else
    if(opname == "S2")
        {
        Op(Up,UpP) = 0.75; 
        Op(Dn,DnP) = 0.75; 
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


inline IQIndexVal SpinHalf::
Up(int i) const
    {
    return getSi(i)(1);
    }

inline IQIndexVal SpinHalf::
Dn(int i) const
    {
    return getSi(i)(2);
    }

inline IQIndexVal SpinHalf::
UpP(int i) const
    {
    return primed(getSi(i))(1);
    }

inline IQIndexVal SpinHalf::
DnP(int i) const
    {
    return primed(getSi(i))(2);
    }


IQTensor inline SpinHalf::
makeSz(int i) const
    {
    IQTensor Sz(conj(si(i)),siP(i));
    Sz(Up(i),UpP(i)) = +0.5;
    Sz(Dn(i),DnP(i)) = -0.5;
    return Sz;
    }

IQTensor inline SpinHalf::
makeSx(int i) const
    {
    IQTensor Sx(conj(si(i)),siP(i));
    Sx(Up(i),DnP(i)) = +0.5;
    Sx(Dn(i),UpP(i)) = +0.5;
    return Sx;
    }

IQTensor inline SpinHalf::
makeISy(int i) const
    {
    IQTensor ISy(conj(si(i)),siP(i));
    ISy(Up(i),DnP(i)) = -0.5;
    ISy(Dn(i),UpP(i)) = +0.5;
    return ISy;
    }

IQTensor inline SpinHalf::
makeSp(int i) const
    {
    IQTensor Sp(conj(si(i)),siP(i));
    Sp(Dn(i),UpP(i)) = 1;
    return Sp;
    }

IQTensor inline SpinHalf::
makeSm(int i) const
    {
    IQTensor Sm(conj(si(i)),siP(i));
    Sm(Up(i),DnP(i)) = 1;
    return Sm;
    }

#endif
