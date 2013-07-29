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

    SpinOne(std::ifstream& s) { doRead(s); }

    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites(const OptSet& opts);
        
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
    if(state == "Up") 
        {
        st = 1;
        }
    else
    if(state == "Z0")
        {
        if(getSi(i).m() == 2)
            Error("Z0 not defined for spin 1/2 site");
        st = 2;
        }
    else
    if(state == "Dn")
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
getOp(int i, const String& opname) const
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
    if(opname == "Sp")
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
    if(opname == "Sm")
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
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

#undef Cout
#undef Endl

#endif
