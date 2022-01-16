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

class BosonSite;

using Boson = BasicSiteSet<BosonSite>;

class BosonSite
    {
    Index s;
    public:

    BosonSite(Index I) : s(I) { }

    BosonSite(Args const& args = Args::global())
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveNb = args.getBool("ConserveNb",conserveQNs);

        auto tags = TagSet("Site,Boson");
        auto n = 1;
        if(args.defined("SiteNumber") )
            {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }

        auto maxOcc = args.getInt("MaxOcc",1);
        if(conserveQNs)
            {
            if(conserveNb)
                {
                auto qints = Index::qnstorage(1+maxOcc);
                for(int n : range(1+maxOcc)) 
                    {
                    qints[n] = QNInt(QN({"Nb",n}),1);
                    }
                s = Index(std::move(qints),tags);
                }
            else
                {
                s = Index(QN(),1+maxOcc,tags);
                }
            }
        else
            {
            if(conserveNb) throw ITError("ConserveNb cannot be true when ConserveQNs=false");
            s = Index(1+maxOcc,tags);
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        auto maxOcc = dim(s)-1;
        for(auto n : range(1+maxOcc))
            {
            if(state == str(n)) return s(1+n);
            }
	if(state == "Emp")return s(1);
        throw ITError("State " + state + " not recognized");
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);
        auto maxOcc = dim(s)-1;

        auto Op = ITensor(dag(s),sP);

        if(opname == "N" || opname == "n")
            {
            for(auto n : range1(1+maxOcc))
                {
                Op.set(s=n,sP=n,n-1);
                }
            }
        else
        if(opname == "A")
            {
            for(auto n : range1(maxOcc))
                {
                Op.set(s=1+n,sP=n,std::sqrt(n));
                }
            }
        else
        if(opname == "Adag")
            {
            for(auto n : range1(maxOcc))
                {
                Op.set(s=n,sP=1+n,std::sqrt(n));
                }
            }
        else
            {
            throw ITError("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

} //namespace itensor

