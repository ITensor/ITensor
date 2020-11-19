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

class SpinHalfSite;

using SpinHalf = BasicSiteSet<SpinHalfSite>;

class SpinHalfSite
    {
    Index s;
    public:

    SpinHalfSite(Index const& I) : s(I) { }

    SpinHalfSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,S=1/2");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        auto conserveqns = args.getBool("ConserveQNs",true);
        auto conserveSz = args.getBool("ConserveSz",conserveqns);
        auto conserveParity = args.getBool("ConserveParity",false);
        if(conserveSz && conserveParity)
            {
            s = Index(QN({"Sz",+1},{"Parity",1,2}),1,
                      QN({"Sz",-1},{"Parity",0,2}),1,Out,ts);
            }
        else if(conserveSz)
            {
            s = Index(QN({"Sz",+1}),1,
                      QN({"Sz",-1}),1,Out,ts);
            }
        else if(conserveParity)
            {
            s = Index(QN({"Parity",1,2}),1,
                      QN({"Parity",0,2}),1,Out,ts);
            }
        else
            {
            s = Index(2,ts);
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "Up") 
            {
            return s(1);
            }
        else 
        if(state == "Dn") 
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
	   Args const& args = Args::global()) const
        {
        auto sP = prime(s);

        auto Up = s(1);
        auto UpP = sP(1);
        auto Dn = s(2);
        auto DnP = sP(2);

        auto Op = ITensor(dag(s),sP);

        if(opname == "Sz")
            {
            Op.set(Up,UpP,+0.5);
            Op.set(Dn,DnP,-0.5);
            }
        else
        if(opname == "Sx")
            {
            //if(not hasQNs(s))
            //    {
                Op.set(Up,DnP,+0.5);
                Op.set(Dn,UpP,+0.5);
            //    }
            //else
            //    {
            //    throw ITError("Operator " + opname + " does not have a well defined QN flux");
            //    }
            }
        else
        if(opname == "ISy")
            {
            //if(not hasQNs(s))
            //    {
                Op.set(Up,DnP,-0.5);
                Op.set(Dn,UpP,+0.5);
            //    }
            //else
            //    {
            //    throw ITError("Operator " + opname + " does not have a well defined QN flux");
            //    }
            }
        else
        if(opname == "Sy")
            {
            //if(not hasQNs(s))
            //    {
                Op.set(Up,DnP,+0.5*Cplx_i);
                Op.set(Dn,UpP,-0.5*Cplx_i);
            //    }
            //else
            //    {
            //    throw ITError("Operator " + opname + " does not have a well defined QN flux");
            //    }
            }
        else
        if(opname == "Sp" || opname == "S+")
            {
            Op.set(Dn,UpP,1);
            }
        else
        if(opname == "Sm" || opname == "S-")
            {
            Op.set(Up,DnP,1);
            }
        else
        if(opname == "projUp")
            {
            Op.set(Up,UpP,1);
            }
        else
        if(opname == "projDn")
            {
            Op.set(Dn,DnP,1);
            }
        else
        if(opname == "S2")
            {
            Op.set(Up,UpP,0.75);
            Op.set(Dn,DnP,0.75);
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

    SpinHalfSite(int n, Args const& args = Args::global())
        {
        *this = SpinHalfSite({args,"SiteNumber=",n});
        }

    };

} //namespace itensor
