//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_Z3_H
#define __ITENSOR_Z3_H
#include "itensor/mps/siteset.h"

namespace itensor {

class Z3 : public SiteSet
    {
    int N_;
    std::vector<IQIndex> site_;
    public:

    Z3();

    Z3(int N);

    Complex static
    Omega()
        {
        static Complex w(cos(2.*Pi/3.),sin(2.*Pi/3.));
        return w;
        }

    private:

    int
    getN() const;

    IQIndex const&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, String const& state) const;

    virtual IQTensor
    getOp(int i, String const& opname, Args const& args) const;

    void
    doRead(std::istream& s);

    void
    doWrite(std::ostream& s) const;

    void
    constructSites();
        
    };

inline Z3::
Z3()
    : N_(-1)
    { }

inline Z3::
Z3(int N)
    : 
    N_(N),
    site_(N_+1)
    { 
    constructSites();
    }

void inline Z3::
constructSites()
    {
    for(int i = 1; i <= N_; ++i)
        {
        site_.at(i) = IQIndex(nameint("Z3 site=",i),
        Index(nameint("0|site",i),1,Site),clock(0,3),
        Index(nameint("1|site",i),1,Site),clock(1,3),
        Index(nameint("2|site",i),1,Site),clock(2,3));
        }
    }

void inline Z3::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

void inline Z3::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

int inline Z3::
getN() const
    { return N_; }

inline 
IQIndex const& Z3::
getSi(int i) const
    { return site_.at(i); }

inline IQIndexVal Z3::
getState(int i, String const& state) const
    {
    int st = -1;
    if(state == "0") 
        {
        st = 1;
        }
    else
    if(state == "1")
        {
        st = 2;
        }
    else
    if(state == "2")
        {
        st = 3;
        }
    else
        {
        Error("State " + state + " not recognized");
        }
    return getSi(i)(st);
    }

inline IQTensor Z3::
getOp(int i, String const& opname, Args const& args) const
    {
    auto s = si(i);
    auto sP = prime(s);

    IQIndexVal Zer(s(1)),
               ZerP(sP(1)),
               One(s(2)),
               OneP(sP(2)),
               Two(s(3)),
               TwoP(sP(3));

    IQTensor Op(dag(s),sP);

    if(opname == "N")
        {
        Op.set(One,OneP,1);
        Op.set(Two,TwoP,2);
        }
    else
    if(opname == "Sig")
        {
        Op.set(Zer,TwoP,1);
        Op.set(One,ZerP,1);
        Op.set(Two,OneP,1);
        }
    else
    if(opname == "SigDag")
        {
        Op.set(Two,ZerP,1);
        Op.set(Zer,OneP,1);
        Op.set(One,TwoP,1);
        }
    else
    if(opname == "Tau")
        {
        Op.set(Zer,ZerP,1);
        Op.set(One,OneP,cos(2.*Pi/3.));
        Op.set(Two,TwoP,cos(4.*Pi/3.));

        IQTensor TauI(s,sP);
        TauI(One,OneP,sin(2.*Pi/3.));
        TauI(Two,TwoP,sin(4.*Pi/3.));

        Op += TauI*Complex_i;
        }
    else
    if(opname == "TauDag")
        {
        Op.set(Zer,ZerP,1);
        Op.set(One,OneP,cos(2.*Pi/3.));
        Op.set(Two,TwoP,cos(4.*Pi/3.));

        IQTensor TauI(s,sP);
        TauI(One,OneP,-sin(2.*Pi/3.));
        TauI(Two,TwoP,-sin(4.*Pi/3.));

        Op += TauI*Complex_i;
        }
    else
    if(opname == "Proj0")
        {
        Op.set(Zer,ZerP,1);
        }
    else
    if(opname == "Proj1")
        {
        Op.set(One,OneP,1);
        }
    else
    if(opname == "Proj2")
        {
        Op.set(Two,TwoP,1);
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

} //namespace itensor

#endif
