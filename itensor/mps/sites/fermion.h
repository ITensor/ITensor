//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_FERMION_H
#define __ITENSOR_FERMION_H
#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class FermionSite;

using Fermion = BasicSiteSet<FermionSite>;

class FermionSite
    {
    Index s;
    public:

    FermionSite() { }

    FermionSite(Index I) : s(I) { }

    FermionSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Fermion");
        auto n = 1;
        if( args.defined("SiteNumber") )
          {
          n = args.getInt("SiteNumber");
          ts.addTags("n="+str(n));
          }
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserve_Nf = args.getBool("ConserveNf",conserveQNs);
        auto oddevenupdown = args.getBool("OddEvenUpDown",false);
        if(not oddevenupdown) //usual case
            {
            auto q_occ = QN({"Nf",1});
            if(not conserve_Nf) q_occ = QN({"Pf",1,-2});
            s = Index(QN(),1,
                      q_occ,1,Out,ts);
            }
        else
            {
            QN q_occ;
            if(n%2==1) q_occ = QN({"Sz",+1},{"Nf",1,-1});
            else       q_occ = QN({"Sz",-1},{"Nf",1,-1});
            s = Index(QN(),1,
                      q_occ,1,Out,ts);
            }
        }

    Index
    index() const { return s; }

    IndexVal
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
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Emp  = s(1);
        auto EmpP = sP(1);
        auto Occ  = s(2);
        auto OccP = sP(2);
         
        auto Op = ITensor(dag(s),sP);

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

    //
    // Deprecated, for backwards compatibility
    //

    FermionSite(int n, Args const& args = Args::global())
        {
        *this = FermionSite({args,"SiteNumber=",n});
        }

    };

//
// Deprecated, for backwards compatability
//

using SpinlessSite = FermionSite;

using Spinless = BasicSiteSet<SpinlessSite>;

} //namespace itensor

#endif
