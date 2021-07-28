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

class Z3Site;

using Z3 = BasicSiteSet<Z3Site>;

class Z3Site
    {
    Index s;
    public:

    Z3Site(Index I) : s(I) { }

    Z3Site(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Z3");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        if(args.getBool("ConserveQNs",true))
          {
          s = Index(QN({"T",0,3}),1,
                    QN({"T",1,3}),1,
                    QN({"T",2,3}),1,Out,ts);
          }
        else
          {
          s = Index(3,ts);
          }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "0") { return s(1); }
        else
        if(state == "1") { return s(2); }
        else
        if(state == "2") { return s(3); }
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

        auto Zer = s(1);
        auto ZerP = sP(1);
        auto One = s(2);
        auto OneP = sP(2);
        auto Two = s(3);
        auto TwoP = sP(3);

        auto Op = ITensor(dag(s),sP);

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
            Op.set(One,OneP,cos(2.*Pi/3.)+sin(2.*Pi/3.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/3.)+sin(4.*Pi/3.)*1_i);
            }
        else
        if(opname == "TauDag")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/3.)-sin(2.*Pi/3.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/3.)-sin(4.*Pi/3.)*1_i);
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
            throw ITError("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    Z3Site(int n, Args const& args = Args::global())
        {
        *this = Z3Site({args,"SiteNumber=",n});
        }

    };

} //namespace itensor

