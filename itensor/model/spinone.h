//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINONE_H
#define __ITENSOR_SPINONE_H
#include "../model.h"

#define Cout std::cout
#define Endl std::endl

class SpinOne : public Model
    {
    public:

    SpinOne();

    SpinOne(int N, const OptSet& opts = Global::opts());

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
    constructSites(const OptSet& opts);
        
    //Older interface for backwards compatibility
    public:
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
    IQTensor
    makeSz2(int i) const;
    IQTensor
    makeSx2(int i) const;
    IQTensor
    makeSy2(int i) const;

    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;

    };

inline SpinOne::
SpinOne()
    : N_(-1)
    { }

inline SpinOne::
SpinOne(int N, const OptSet& opts)
    : N_(N),
      site_(N_+1)
    { 
    constructSites(opts);
    }

inline void SpinOne::
constructSites(const OptSet& opts)
    {
    for(int j = 1; j <= N_; ++j)
        {
        if((opts.getBool("SHalfEdge",false) && (j==1 || j == N_))
           || (opts.getBool("SHalfLeftEdge",false) && j==1))
            {
            if(opts.getBool("Verbose",false))
                {
                Cout << "Placing a S=1/2 at site " << j << Endl;
                }

            site_.at(j) = IQIndex(nameint("S=1/2 site=",j),
                Index(nameint("Up for site",j),1,Site),QN(+1,0),
                Index(nameint("Dn for site",j),1,Site),QN(-1,0));
            }
        else
            {
            site_.at(j) = IQIndex(nameint("S=1 site=",j),
                Index(nameint("Up for site",j),1,Site),QN(+2,0),
                Index(nameint("Z0 for site",j),1,Site),QN( 0,0),
                Index(nameint("Dn for site",j),1,Site),QN(-2,0));
            }
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
getN() const
    { return N_; }

inline const IQIndex& SpinOne::
getSi(int i) const
    { return site_.at(i); }

inline IQIndexVal SpinOne::
getState(int i, const String& state) const
    {
    int st = -1;
    if(state == "Up" || state == "+") 
        {
        st = 1;
        }
    else
    if(state == "Z0" || state == "0")
        {
        if(getSi(i).m() == 2)
            Error("Z0 not defined for spin 1/2 site");
        st = 2;
        }
    else
    if(state == "Dn" || state == "-")
        {
        st = getSi(i).m();
        }
    else
        {
        Error("State " + state + " not recognized");
        }
    return getSi(i)(st);
    }

inline IQTensor SpinOne::
getOp(int i, const String& opname, const OptSet& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = primed(s);

    IQIndexVal Up(s(1)),
               UpP(sP(1)),
               Dn(s(s.m())),
               DnP(sP(s.m())),
               Z0(s(2)),
               Z0P(sP(2));

    IQTensor Op(conj(s),sP);

    if(opname == "Sz")
        {
        if(s.m() == 2)
            {
            Op(Up,UpP) = +0.5;
            Op(Dn,DnP) = -0.5;
            }
        else
            {
            Op(Up,UpP) = +1.;
            Op(Dn,DnP) = -1.;
            }
        }
    else
    if(opname == "Sx")
        {
        if(s.m() == 2)
            {
            Op(Up,DnP) = +0.5;
            Op(Dn,UpP) = +0.5;
            }
        else
            {
            Op(Up,Z0P) = ISqrt2; 
            Op(Z0,UpP) = ISqrt2;
            Op(Z0,DnP) = ISqrt2; 
            Op(Dn,Z0P) = ISqrt2;
            }
        }
    else
    if(opname == "ISy")
        {
        if(s.m() == 2)
            {
            Op(Up,DnP) = -0.5;
            Op(Dn,UpP) = +0.5;
            }
        else
            {
            Op(Up,Z0P) = +ISqrt2; 
            Op(Z0,UpP) = -ISqrt2;
            Op(Z0,DnP) = +ISqrt2; 
            Op(Dn,Z0P) = -ISqrt2;
            }
        }
    else
    if(opname == "Sy")
        {
        if(s.m() == 2)
            {
            Op(Up,DnP) = +0.5;
            Op(Dn,UpP) = -0.5;
            }
        else
            {
            Op(Up,Z0P) = -ISqrt2; 
            Op(Z0,UpP) = +ISqrt2;
            Op(Z0,DnP) = -ISqrt2; 
            Op(Dn,Z0P) = +ISqrt2;
            }
        Op *= Complex_i;
        }
    else
    if(opname == "Sp" || opname == "S+")
        {
        if(s.m() == 2)
            {
            Op(Dn,UpP) = 1;
            }
        else
            {
            Op(Dn,Z0P) = Sqrt2; 
            Op(Z0,UpP) = Sqrt2;
            }
        }
    else
    if(opname == "Sm" || opname == "S-")
        {
        if(s.m() == 2)
            {
            Op(Up,DnP) = 1;
            }
        else
            {
            Op(Up,Z0P) = Sqrt2; 
            Op(Z0,DnP) = Sqrt2;
            }
        }
    else
    if(opname == "Sz2")
        {
        if(s.m() == 2) Error("Sz^2 only non-trivial for S=1 sites");
        Op(Up,UpP) = 1; 
        Op(Dn,DnP) = 1;
        }
    else
    if(opname == "Sx2")
        {
        if(s.m() == 2) Error("Sx^2 only non-trivial for S=1 sites");
        Op(Up,UpP) = 0.5; 
        Op(Up,DnP) = 0.5;
        Op(Z0,Z0P) = 1;
        Op(Dn,DnP) = 0.5; 
        Op(Dn,UpP) = 0.5;
        }
    else
    if(opname == "Sy2")
        {
        if(s.m() == 2) Error("Sy^2 only non-trivial for S=1 sites");
        Op(Up,UpP) = 0.5; 
        Op(Up,DnP) = -0.5;
        Op(Z0,Z0P) = 1;
        Op(Dn,DnP) = 0.5; 
        Op(Dn,UpP) = -0.5;
        }
    else
    if(opname == "projUp")
        {
        Op(Up,UpP) = 1; 
        }
    else
    if(opname == "projZ0")
        {
        if(s.m() == 2) Error("Can only form projZ0 for S=1 sites");
        Op(Z0,Z0P) = 1; 
        }
    else
    if(opname == "projDn")
        {
        Op(Dn,DnP) = 1; 
        }
    else
    if(opname == "XUp")
        {
        //m = +1 state along x axis
        Op = IQTensor(s);
        Op(Up) = 0.5;
        Op(Z0) = ISqrt2;
        Op(Dn) = 0.5;
        }
    else
    if(opname == "XZ0")
        {
        //m = 0 state along x axis
        Op = IQTensor(s);
        Op(Up) = ISqrt2;
        Op(Dn) = -ISqrt2;
        }
    else
    if(opname == "XDn")
        {
        //m = -1 state along x axis
        Op = IQTensor(s);
        Op(Up) = 0.5;
        Op(Z0) = -ISqrt2;
        Op(Dn) = 0.5;
        }
    else
    if(opname == "S2")
        {
        const Real ssp1 = (s.m()==2 ? 0.75 : 2.);
        Op(Up,UpP) = ssp1; 
        Op(Dn,DnP) = ssp1;
        if(s.m() > 2)
            Op(Z0,Z0P) = ssp1;
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
    return primed(getSi(i))(1);
    }

inline IQIndexVal SpinOne::
Z0P(int i) const
    {
    IQIndex ind = primed(getSi(i));
    if(ind.m() == 2)
        Error("Z0 not defined for spin 1/2 site");
    return ind(2);
    }

inline IQIndexVal SpinOne::
DnP(int i) const
    {
    IQIndex ind = primed(getSi(i));
    return ind(ind.m());
    }

IQTensor inline SpinOne::
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

IQTensor inline SpinOne::
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

IQTensor inline SpinOne::
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

IQTensor inline SpinOne::
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

IQTensor inline SpinOne::
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


IQTensor inline SpinOne::
makeSz2(int i) const
    {
    if(si(i).m() == 2) Error("Sz^2 only non-trivial for S=1 sites");
    IQTensor Sz2(conj(si(i)),siP(i));
    Sz2(Up(i),UpP(i)) = 1; 
    Sz2(Dn(i),DnP(i)) = 1;
    return Sz2;
    }

IQTensor inline SpinOne::
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

IQTensor inline SpinOne::
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

#undef Cout
#undef Endl

#endif
