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
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/util/str.h"

// Code written by Samuel Gozel

namespace itensor {

class SpinTwo : public SiteSet
    {
    public:

    SpinTwo() { }

    SpinTwo(int N, 
        Args const& args = Args::global());

    void
        read(std::istream& s);

    };


class SpinTwoSite
	{
  Index s;
	public:

    SpinTwoSite(Index I) : s(I) { }

    SpinTwoSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,S=2");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);
        if(conserveSz)
            {
            s = Index{QN({"Sz",+4}),1,
                      QN({"Sz",+2}),1,
                      QN({"Sz",0}),1,
                      QN({"Sz",-2}),1,
                      QN({"Sz",-4}),1,Out,ts};
            }
        else
            {
            s = Index{5,ts};
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
		{
        if (state == "Up" || state == "4") 
            {
            return s(1);
            }
        else if (state == "Upi" || state == "2")
            {
            return s(2);
            }
        else if (state == "Z0" || state == "0")
            {
            return s(3);
            }
        else if (state == "Dni" || state == "-2")
            {
            return s(4);
            }
        else if (state == "Dn" || state == "-4")
            {
            return s(5);
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
        const Real val1 = std::sqrt(6.0)/2.0;
        const Real val2 = std::sqrt(6.0);

        auto sP = prime(s);

        auto Up  = s(1);
        auto UpP = sP(1);
        auto Upi = s(2);
        auto UpiP = sP(2);
        auto Z0  = s(3);
        auto Z0P = sP(3);
        auto Dni = s(4);
        auto DniP = sP(4);
        auto Dn  = s(5);
        auto DnP = sP(5);

        auto Op = ITensor(dag(s),sP);

        if (opname == "Sz")
            {
            Op.set(Up,UpP,+2.0);
            Op.set(Upi,UpiP,+1.0);
            Op.set(Dni,DniP,-1.0);
            Op.set(Dn,DnP,-2.0);
            }
        else if (opname == "Sx")
            {
            Op.set(Up,UpiP,1.0);
            Op.set(Upi,UpP,1.0);
            Op.set(Upi,Z0P,val1); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Z0,UpiP,val1);
            Op.set(Z0,DniP,val1);
            Op.set(Dni,Z0P,val1);
            Op.set(Dni,DnP,1.0);
            Op.set(Dn,DniP,1.0);
            }
        else if (opname == "ISy") // defined as i*Sy
            {
            Op.set(Up,UpiP,-1.0);
            Op.set(Upi,UpP,1.0);
            Op.set(Upi,Z0P,-val1);
            Op.set(Z0,UpiP,val1);
            Op.set(Z0,DniP,-val1);
            Op.set(Dni,Z0P,val1);
            Op.set(Dni,DnP,-1.0);
            Op.set(Dn,DniP,1.0);
            }
        else if (opname == "Sy")
            {
            Op.set(Up,UpiP,1.0*Cplx_i);
            Op.set(Upi,UpP,-1.0*Cplx_i);
            Op.set(Upi,Z0P,val1*Cplx_i);
            Op.set(Z0,UpiP,-val1*Cplx_i);
            Op.set(Z0,DniP,val1*Cplx_i);
            Op.set(Dni,Z0P,-val1*Cplx_i);
            Op.set(Dni,DnP,1.0*Cplx_i);
            Op.set(Dn,DniP,-1.0*Cplx_i);
            }
        else if (opname == "Sp" || opname == "S+")
            {
            Op.set(Upi,UpP,2.0);
            Op.set(Z0,UpiP,val2);
            Op.set(Dni,Z0P,val2);
            Op.set(Dn,DniP,2.0);
            }
        else if (opname == "Sm" || opname == "S-")
            {
            Op.set(Up,UpiP,2.0);
            Op.set(Upi,Z0P,val2);
            Op.set(Z0,DniP,val2);
            Op.set(Dni,DnP,2.0);
            }
        else if (opname == "Sz2")
            {
            Op.set(Up,UpP,4);
            Op.set(Upi,UpiP,1);
            Op.set(Dni,DniP,1);
            Op.set(Dn,DnP,4);
            }
        else if (opname == "Sx2")
            {
            Op.set(Up,UpP,1.0);
            Op.set(Up,Z0P,val1);
            Op.set(Upi,UpiP,2.5);
            Op.set(Upi,DniP,1.5);
            Op.set(Z0,UpP,val1);
            Op.set(Z0,Z0P,3.0);
            Op.set(Z0,DnP,val1);
            Op.set(Dni,UpiP,1.5);
            Op.set(Dni,DniP,2.5);
            Op.set(Dn,Z0P,val1);
            Op.set(Dn,DnP,1.0);
            }
        else if (opname == "Sy2")
            {
            Op.set(Up,UpP,1.0);
            Op.set(Up,Z0P,-val1);
            Op.set(Upi,UpiP,2.5);
            Op.set(Upi,DniP,-1.5);
            Op.set(Z0,UpP,-val1);
            Op.set(Z0,Z0P,3.0);
            Op.set(Z0,DnP,-val1);
            Op.set(Dni,UpiP,-1.5);
            Op.set(Dni,DniP,2.5);
            Op.set(Dn,Z0P,-val1);
            Op.set(Dn,DnP,1.0);
            }
        else if (opname == "projUp")
            {
            Op.set(Up,UpP,1);
            }
        else if (opname == "projUpi")
            {
            Op.set(Upi,UpiP,1);
            }
        else if (opname == "projZ0")
            {
            Op.set(Z0,Z0P,1);
            }
        else if (opname == "projDni")
            {
            Op.set(Dni,DniP,1);
            }
        else if (opname == "projDn")
            {
            Op.set(Dn,DnP,1);
            }
        else if (opname == "S2")
            {
            Op.set(Up,UpP,6);
            Op.set(Upi,UpiP,6);
            Op.set(Z0,Z0P,6);
            Op.set(Dni,DniP,6);
            Op.set(Dn,DniP,6);
            }
        else
            {
            throw ITError("Operator " + opname + " name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    SpinTwoSite(int n, Args const& args = Args::global())
        {
        *this = SpinTwoSite({args,"SiteNumber=",n});
        }

    }; //SpinTwoSite

inline SpinTwo::
SpinTwo(int N, 
        Args const& args)
	{
    auto shedge = args.getBool("SHalfEdge",false);
    auto Lshedge = args.getBool("SHalfLeftEdge",false);

    auto sites = SiteStore(N);

    auto start = 1;
    if(shedge || Lshedge)
        {
        if(args.getBool("Verbose",false)) println("Placing a S=1/2 at site 1");
        sites.set(1,SpinHalfSite({args,"SiteNumber=",1}));
        start = 2;
        }

    for(int j = start; j < N; ++j)
        {
        sites.set(j,SpinTwoSite({args,"SiteNumber=",j}));
        }

    if(shedge)
        {
        if(args.getBool("Verbose",false)) println("Placing a S=1/2 at site N=",N);
        sites.set(N,SpinHalfSite({args,"SiteNumber=",N}));
        }
    else
        {
        sites.set(N,SpinTwoSite({args,"SiteNumber=",N}));
        }

    SiteSet::init(std::move(sites));
	}

void inline SpinTwo::
read(std::istream& s)
	{
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        for(int j = 1; j <= N; ++j) 
            {
            auto I = Index{};
            I.read(s);
            if(dim(I) == 5) store.set(j,SpinTwoSite(I));
            else if(dim(I) == 2) store.set(j,SpinHalfSite(I));
            else throw ITError(format("SpinTwo cannot read index of size %d",dim(I)));
            }
        init(std::move(store));
        }
	}

} //namespace itensor

