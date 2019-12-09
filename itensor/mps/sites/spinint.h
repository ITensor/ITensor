//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINI_H
#define __ITENSOR_SPINI_H
#include "itensor/mps/siteset.h"

namespace itensor {

class SpinSSite;

using SpinInt = BasicSiteSet<SpinSSite>;

// class SpinS : public SiteSet
//     {
//     public:
// 
//     SpinS() { }
// 
//     SpinS(int N, 
//             Args const& args = Args::global());
// 
//     void
//     read(std::istream& s);
// 
//     };


class SpinSSite
    {
    Index s;
    //std::vector<std::string> Occn;
    public:

    SpinSSite() { }

    SpinSSite(Index I) : s(I) { }
    
    SpinSSite(Args const& args = Args::global())
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        //auto conserveNb = args.getBool("ConserveNb",conserveQNs);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);

        auto tags = TagSet("Site,SpinInt");
        auto n = 1;
        if(args.defined("SiteNumber") )
            {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }

        auto Smax = args.getInt("Smax",1);
        if(conserveQNs)
            {
            if(conserveSz)
                {
                auto qints = Index::qnstorage(2*Smax+1);
                for(int n : range(2*Smax+1)) 
                    {
                    qints[n] = QNInt(QN({"Sz",static_cast<QNum::qn_t>(2*(n-Smax))}),1);
                    }
                s = Index(std::move(qints),tags);
                }
            else
                {
                s = Index(QN(),2*Smax+1,tags);
                }
            }
        else
            {
            if(conserveSz) Error("ConserveSz cannot be true when ConserveQNs=false");
            s = Index(2*Smax+1,tags);
            }
        }
    
    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        auto Smax = (dim(s)-1)/2;
        for(auto n : range(2*Smax+1))
            {
            if(state == str(n-Smax)) return s(1+n);
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
        
        auto Smax = (dim(s)-1)/2;

        if(opname == "Sz")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -Smax+i-1;
                Op.set(s(i),sP(i),sz);
                }
            }
        else
        if(opname == "S+")
            {
            //Op = mixedIQTensor(s,sP);
            for(int i = 1; i <= dim(s)-1; ++i)
                {
                auto sz = -Smax+i-1;
                Op.set(s(i),sP(i+1),std::sqrt((Smax-sz)*(Smax+sz+1.0)));
                }
            }
        else
        if(opname == "S-")
            {
            for(int i = 2; i <= dim(s); ++i)
                {
                auto sz = -Smax+i-1;
                Op.set(s(i),sP(i-1),std::sqrt((Smax+sz)*(Smax-sz+1.0)));
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
                Op.set(s(i),sP(i-1),1.0);
                }
            }
        else
        if(opname == "Sz2")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -Smax+i-1;
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
