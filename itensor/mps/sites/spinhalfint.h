//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINHI_H
#define __ITENSOR_SPINHI_H
#include "itensor/mps/siteset.h"

namespace itensor {

class SpinShSite;

using SpinHalfInt = BasicSiteSet<SpinShSite>;


class SpinShSite
    {
    Index s;
    public:

    SpinShSite() { }

    SpinShSite(Index I) : s(I) { }
    
    SpinShSite(Args const& args = Args::global())
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);

        auto tags = TagSet("Site,SpinHalfInt");
        auto n = 1;
        if(args.defined("SiteNumber") )
            {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }

        auto DSmax = args.getInt("DSmax",3); //2*Smax
        if(conserveQNs)
            {
            if(conserveSz)
                {
                auto qints = Index::qnstorage(DSmax+1);
                for(int n : range(DSmax+1)) 
                    {
                    qints[n] = QNInt(QN({"Sz",static_cast<QNum::qn_t>(2*n-DSmax)}),1);
                    }
                s = Index(std::move(qints),tags);
                }
            else
                {
                s = Index(QN(),DSmax+1,tags);
                }
            }
        else
            {
            if(conserveSz) Error("ConserveSz cannot be true when ConserveQNs=false");
            s = Index(DSmax+1,tags);
            }
        }


    Index
    index() const { return s; }
    
    IndexVal
    state(std::string const& state)
        {
        auto DSmax = dim(s)-1;
        for(auto n : range(DSmax+1))
            {
            if(state == str(2*n-DSmax)) return s(1+n); //state name is "2Sz"
            }
        Error("State " + state + " not recognized");
        return IndexVal{};
        }


	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Op = ITensor(dag(s),sP);
        
        auto DSmax = dim(s)-1;

        if(opname == "Sz")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i),sz);
                }
            }
        else
        if(opname == "S+")
            {
            for(int i = 1; i <= dim(s)-1; ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i+1),std::sqrt((DSmax/2.0-sz)*(DSmax/2.0+sz+1.0)));
                }
            }
        else
        if(opname == "S-")
            {
            for(int i = 2; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i-1),std::sqrt((DSmax/2.0+sz)*(DSmax/2.0-sz+1.0)));
                }
            }
        else
        if(opname == "U+")
            {
            for(int i = 1; i <= dim(s)-1; ++i)
                {
                Op.set(s(i),sP(i+1),1.0);
                }
            }
        else
        if(opname == "U-")
            {
            for(int i = 2; i <= dim(s); ++i)
                {
                //auto sz = -(DSpinS-1)/2+i-1;
                Op.set(s(i),sP(i-1),1.0);
                }
            }
        else
        if(opname == "Sz2")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i),pow(sz, 2.0));
                }
            }
        else
            {
            Error("Operator " + opname + " name not recognized");
            }

        return Op;
        }
    };
}
#endif
