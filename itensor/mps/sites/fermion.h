//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#pragma once

#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class FermionSite;

using Fermion = BasicSiteSet<FermionSite>;

class FermionSite
    {
    Index s;
    public:

    FermionSite(Index I) : s(I) { }

    FermionSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Fermion");
        auto n = 1;
        if(args.defined("SiteNumber"))
          {
          n = args.getInt("SiteNumber");
          ts.addTags("n="+str(n));
          }
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserve_Nf = args.getBool("ConserveNf",conserveQNs);
        auto oddevenupdown = args.getBool("OddEvenUpDown",false);
        if(not conserveQNs)
            {
            s = Index(2,ts);
            }
        else if(not oddevenupdown) //usual case
            {
            if(conserve_Nf) //usual case
                {
                s = Index(QN({"Nf",0,-1}),1,
                          QN({"Nf",1,-1}),1,Out,ts);
                }
            else
                {
                s = Index(QN({"Pf",0,-2}),1,
                          QN({"Pf",1,-2}),1,Out,ts);
                }
            }
        else
            {
            auto q_emp = QN({"Sz",0},{"Nf",0,-1});
            QN q_occ;
            if(n%2==1) q_occ = QN({"Sz",+1},{"Nf",1,-1});
            else       q_occ = QN({"Sz",-1},{"Nf",1,-1});
            s = Index(q_emp,1,
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
            throw ITError("State " + state + " not recognized");
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
            throw ITError("Operator \"" + opname + "\" name not recognized");
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

