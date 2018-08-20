//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINLESS_H
#define __ITENSOR_SPINLESS_H
#include "itensor/mps/siteset.h"

namespace itensor {

class SpinlessSite;

using Spinless = BasicSiteSet<SpinlessSite>;

class SpinlessSite
    {
    IQIndex s;
    public:

    SpinlessSite() { }

    SpinlessSite(IQIndex I) : s(I) { }

    SpinlessSite(int n, Args const& args = Args::global())
        {
        auto conserve_Nf = args.getBool("ConserveNf",true);
        auto oddevenupdown = args.getBool("OddEvenUpDown",false);

        if(!oddevenupdown) //usual case
            {
            auto q_occ = QN("Nf=",1);
            if(not conserve_Nf) q_occ = QN("Pf=",1);
            s = IQIndex{nameint("Spinless ",n),
                Index(nameint("Emp ",n),1,Site),QN(),
                Index(nameint("Occ ",n),1,Site),q_occ};
            }
        else
            {
            QN q_occ;
            if(n%2==1) q_occ = QN("Sz",+1,"Nf=",1);
            else       q_occ = QN("Sz",-1,"Nf=",1);
            s = IQIndex{nameint("Spinless ",n),
                Index(nameint("Emp ",n),1,Site),QN(),
                Index(nameint("Occ ",n),1,Site),q_occ};
            }
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "Emp" || state == "0") 
            {
            return s(1);
            }
        else 
        if(state == "Occ" || state == "1") 
            {
            return s(2);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

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
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

} //namespace itensor

#endif
