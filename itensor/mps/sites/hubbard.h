//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HUBBARD_H
#define __ITENSOR_HUBBARD_H
#include "../siteset.h"

namespace itensor {

class Hubbard : public SiteSet
    {
    public:

    Hubbard();

    Hubbard(int N, 
            const Args& opts = Global::opts());

    bool
    conserveNf() const { return conserveNf_; }

    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const Args& opts = Global::opts()) const;

    DefaultOpsT
    getDefaultOps(const Args& opts) const;

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
Hubbard(int N, const Args& opts)
    : N_(N),
      site_(N_+1)
    { 
    conserveNf_ = opts.getBool("ConserveNf",true);
    conserveSz_ = opts.getBool("ConserveSz",true);
    constructSites();
    }

void inline Hubbard::
constructSites()
    {
    const int One = (conserveNf_ ? 1 : 0),
              Two = 2*One,
              Up = (conserveSz_ ? +1 : 0),
              Dn = -Up;
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("Hubbard site=",j),
            Index(nameint("Emp for site ",j),1,Site),  QN( 0,0,0),
            Index(nameint("Up for site ",j),1,Site),   QN(Up,One,1),
            Index(nameint("Dn for site ",j),1,Site),   QN(Dn,One,1),
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
    return IQIndexVal();
    }


inline IQTensor Hubbard::
getOp(int i, const String& opname, const Args& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = prime(s);

    IQIndexVal Em(s(1)),
               EmP(sP(1)),
               Up(s(2)),
               UpP(sP(2)),
               Dn(s(3)),
               DnP(sP(3)),
               UD(s(4)),
               UDP(sP(4));

    IQTensor Op(dag(s),sP);

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
    if(opname == "S2")
        {
        //S dot S on-site
        Op(Up,UpP) = 0.75; 
        Op(Dn,DnP) = 0.75;
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

Hubbard::DefaultOpsT inline Hubbard::
getDefaultOps(const Args& opts) const
    {
    static const std::vector<String> dops_(initDefaultOps());
    return dops_;
    }

}; //namespace itensor

#endif
