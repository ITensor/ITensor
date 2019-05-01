//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BOSON_H
#define __ITENSOR_BOSON_H
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
        auto tags = TagSet("Site,Boson");
        auto n = 1;
        if(args.defined("SiteNumber") )
            {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }

        auto maxOcc = args.getInt("MaxOcc",1);
        if(args.getBool("ConserveQNs",true))
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
        Error("State " + state + " not recognized");
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
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

} //namespace itensor

#endif
