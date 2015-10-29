//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_Z3_H
#define __ITENSOR_Z3_H
#include "itensor/mps/sites/siteset.h"

namespace itensor {

class Z3 : public SiteSet
    {
    public:

    Z3();

    Z3(int N);


    Complex static
    Omega()
        {
        static Complex w(cos(2.*Pi/3.),sin(2.*Pi/3.));
        return w;
        }


    //Operators

    private:

    int
    getN() const;

    const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const Args& opts) const;

    void
    doRead(std::istream& s);

    void
    doWrite(std::ostream& s) const;

    void
    constructSites();
        
    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;

    };

inline Z3::
Z3()
    : N_(-1)
    { 
    QN::Nmax() = 3;
    }

inline Z3::
Z3(int N)
    : 
    N_(N),
    site_(N_+1)
    { 
    QN::Nmax() = 3;
    constructSites();
    }

void inline Z3::
constructSites()
    {
    QN::Nmax() = 3;
    for(int i = 1; i <= N_; ++i)
        {
        site_.at(i) = IQIndex(nameint("Z3 site=",i),
        Index(nameint("0|site",i),1,Site),QN(0,0,0),
        Index(nameint("1|site",i),1,Site),QN(0,0,1),
        Index(nameint("2|site",i),1,Site),QN(0,0,2));
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
const IQIndex& Z3::
getSi(int i) const
    { return site_.at(i); }

inline IQIndexVal Z3::
getState(int i, const String& state) const
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
getOp(int i, const String& opname, const Args& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = prime(s);

    IQIndexVal Zer(s(1)),
               ZerP(sP(1)),
               One(s(2)),
               OneP(sP(2)),
               Two(s(3)),
               TwoP(sP(3));

    IQTensor Op(dag(s),sP);

    if(opname == "N")
        {
        Op(One,OneP) = 1;
        Op(Two,TwoP) = 2;
        }
    else
    if(opname == "Sig")
        {
        Op(Zer,TwoP) = 1;
        Op(One,ZerP) = 1;
        Op(Two,OneP) = 1;
        }
    else
    if(opname == "SigDag")
        {
        Op(Two,ZerP) = 1;
        Op(Zer,OneP) = 1;
        Op(One,TwoP) = 1;
        }
    else
    if(opname == "Tau")
        {
        Op(Zer,ZerP) = 1;
        Op(One,OneP) = cos(2.*Pi/3.);
        Op(Two,TwoP) = cos(4.*Pi/3.);

        IQTensor TauI(s,sP);
        TauI(One,OneP) = sin(2.*Pi/3.);
        TauI(Two,TwoP) = sin(4.*Pi/3.);

        Op += TauI*Complex_i;
        }
    else
    if(opname == "TauDag")
        {
        Op(Zer,ZerP) = 1;
        Op(One,OneP) = cos(2.*Pi/3.);
        Op(Two,TwoP) = cos(4.*Pi/3.);

        IQTensor TauI(s,sP);
        TauI(One,OneP) = -sin(2.*Pi/3.);
        TauI(Two,TwoP) = -sin(4.*Pi/3.);

        Op += TauI*Complex_i;
        }
    else
    if(opname == "Proj0")
        {
        Op(Zer,ZerP) = 1;
        }
    else
    if(opname == "Proj1")
        {
        Op(One,OneP) = 1;
        }
    else
    if(opname == "Proj2")
        {
        Op(Two,TwoP) = 1;
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

}; //namespace itensor

#endif
