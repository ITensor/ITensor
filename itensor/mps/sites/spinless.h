//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINLESS_H
#define __ITENSOR_SPINLESS_H
#include "itensor/mps/siteset.h"

namespace itensor {

class Spinless : public SiteSet
    {
    public:

    Spinless();

    Spinless(int N, Args const& args = Args::global());

    private:

    virtual IQIndexVal
    getState(int i, String const& state) const;

    virtual IQTensor
    getOp(int i, String const& opname, Args const& args = Args::global()) const;

    void
    constructSites();

    void
    doRead(std::istream& s);

    void
    doWrite(std::ostream& s) const;

    //Data members -----------------

    bool conserve_Nf_;

    };

inline Spinless::
Spinless()
  : conserve_Nf_(true)
    { 
    }

inline Spinless::
Spinless(int N, Args const& args)
  : SiteSet(N)
    { 
    conserve_Nf_ = args.getBool("ConserveNf",true);
    constructSites();
    }

void inline Spinless::
constructSites()
    {
    auto q_occ = QN("Nf=",1);
    if(not conserve_Nf_)
        {
        q_occ = QN("Pf=",1);
        }
    for(int i = 1; i <= N(); ++i)
        {
        set(i,{nameint("Spinless site=",i),
        Index(nameint("Emp for site",i),1,Site),QN(),
        Index(nameint("Occ for site",i),1,Site),q_occ});
        }
    }

void inline Spinless::
doRead(std::istream& s)
    {
    s.read((char*) &conserve_Nf_,sizeof(conserve_Nf_));
    }

void inline Spinless::
doWrite(std::ostream& s) const
    {
    s.write((char*) &conserve_Nf_,sizeof(conserve_Nf_));
    }

inline IQIndexVal Spinless::
getState(int i, String const& state) const
    {
    if(state == "Emp" || state == "0") 
        {
        return si(i)(1);
        }
    else 
    if(state == "Occ" || state == "1") 
        {
        return si(i)(2);
        }
    else
        {
        Error("State " + state + " not recognized");
        }
    return IQIndexVal{};
    }

inline IQTensor Spinless::
getOp(int i, String const& opname, Args const& args) const
    {
    auto s  = si(i);
    auto sP = prime(si(i));

    auto Emp  = s(1);
    auto EmpP = sP(1);
    auto Occ  = s(2);
    auto OccP = sP(2);
     
    auto Op = IQTensor(dag(s),sP);

    if(opname == "N" || opname == "n")
        {
        Op.set(Occ,OccP,1);
        }
    else
    if(opname == "C")
        {
        Op.set(Occ,EmpP,1);
        }
    else
    if(opname == "Cdag")
        {
        Op.set(Emp,OccP,1);
        }
    else
    if(opname == "A")
        {
        Op.set(Occ,EmpP,1);
        }
    else
    if(opname == "Adag")
        {
        Op.set(Emp,OccP,1);
        }
    else
    if(opname == "F" || opname == "FermiPhase")
        {
        Op.set(Emp,EmpP,1);
        Op.set(Occ,OccP,-1);
        }
    else
    if(opname == "projEmp")
        {
        Op.set(Emp,EmpP,1);
        }
    else
    if(opname == "projOcc")
        {
        Op.set(Occ,OccP,1); 
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

} //namespace itensor

#endif
