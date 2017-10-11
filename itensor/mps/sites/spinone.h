//
// Distributed under the ITensor Library License, Version 1.2
// (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINONE_H
#define __ITENSOR_SPINONE_H
#include "itensor/mps/siteset.h"

namespace itensor {

class SpinOne : public SiteSet
    {
    public:

    SpinOne() { }

    SpinOne(int N, 
            Args const& args = Args::global());

    SpinOne(std::vector<IQIndex> const& inds);

    void
    read(std::istream& s);

    };


class SpinOneSite
    {
    IQIndex s;
    public:

    SpinOneSite() { }

    SpinOneSite(IQIndex I) : s(I) { }

    SpinOneSite(int n, Args const& args = Args::global())
        {
        s = IQIndex{nameint("S=1 site=",n),
            Index(nameint("Up:site",n),1,Site),QN("Sz=",+2),
            Index(nameint("Z0:site",n),1,Site),QN("Sz=", 0),
            Index(nameint("Dn:site",n),1,Site),QN("Sz=",-2)};
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "Up" || state == "+") 
            {
            return s(1);
            }
        else
        if(state == "Z0" || state == "0")
            {
            return s(2);
            }
        else
        if(state == "Dn" || state == "-")
            {
            return s(3);
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

        auto Up  = s(1);
        auto UpP = sP(1);
        auto Z0  = s(2);
        auto Z0P = sP(2);
        auto Dn  = s(s.m());
        auto DnP = sP(s.m());

        auto Op = IQTensor(dag(s),sP);

        if(opname == "Sz")
            {
            Op.set(Up,UpP,+1.0);
            Op.set(Dn,DnP,-1.0);
            }
        else
        if(opname == "Sx")
            {
            //mixedIQTensor call needed here
            //because as an IQTensor, Op would
            //not have a well defined QN flux
            Op = mixedIQTensor(dag(s),sP);
            Op.set(Up,Z0P,ISqrt2); 
            Op.set(Z0,UpP,ISqrt2);
            Op.set(Z0,DnP,ISqrt2); 
            Op.set(Dn,Z0P,ISqrt2);
            }
        else
        if(opname == "ISy")
            {
            //mixedIQTensor call needed here
            //because as an IQTensor, Op would
            //not have a well defined QN flux
            Op = mixedIQTensor(dag(s),sP);
            Op.set(Up,Z0P,-ISqrt2); 
            Op.set(Z0,UpP,+ISqrt2);
            Op.set(Z0,DnP,-ISqrt2); 
            Op.set(Dn,Z0P,+ISqrt2);
            }
        else
        if(opname == "Sy")
            {
            //mixedIQTensor call needed here
            //because as an IQTensor, Op would
            //not have a well defined QN flux
            Op = mixedIQTensor(dag(s),sP);
            Op.set(Up,Z0P,+ISqrt2*1_i); 
            Op.set(Z0,UpP,-ISqrt2*1_i);
            Op.set(Z0,DnP,+ISqrt2*1_i); 
            Op.set(Dn,Z0P,-ISqrt2*1_i);
            }
        else
        if(opname == "Sp" || opname == "S+")
            {
            Op.set(Dn,Z0P,Sqrt2);  
            Op.set(Z0,UpP,Sqrt2);
            }
        else
        if(opname == "Sm" || opname == "S-")
            {
            Op.set(Up,Z0P,Sqrt2);
            Op.set(Z0,DnP,Sqrt2);
            }
        else
        if(opname == "Sz2")
            {
            Op.set(Up,UpP,1); 
            Op.set(Dn,DnP,1);
            }
        else
        if(opname == "Sx2")
            {
            Op = mixedIQTensor(dag(s),sP);
            Op.set(Up,UpP,0.5); 
            Op.set(Up,DnP,0.5);
            Op.set(Z0,Z0P,1.0);
            Op.set(Dn,DnP,0.5); 
            Op.set(Dn,UpP,0.5);
            }
        else
        if(opname == "Sy2")
            {
            Op = mixedIQTensor(dag(s),sP);
            Op.set(Up,UpP,+0.5); 
            Op.set(Up,DnP,-0.5);
            Op.set(Z0,Z0P,1);
            Op.set(Dn,DnP,+0.5); 
            Op.set(Dn,UpP,-0.5);
            }
        else
        if(opname == "projUp")
            {
            Op.set(Up,UpP,1);
            }
        else
        if(opname == "projZ0")
            {
            Op.set(Z0,Z0P,1);
            }
        else
        if(opname == "projDn")
            {
            Op.set(Dn,DnP,1);
            }
        else
        if(opname == "XUp")
            {
            //m = +1 state along x axis
            Op = mixedIQTensor(dag(s),sP);
            Op = IQTensor(s);
            Op.set(Up,0.5);
            Op.set(Z0,ISqrt2);
            Op.set(Dn,0.5);
            }
        else
        if(opname == "XZ0")
            {
            //m = 0 state along x axis
            Op = mixedIQTensor(dag(s),sP);
            Op = IQTensor(s);
            Op.set(Up,+ISqrt2);
            Op.set(Dn,-ISqrt2);
            }
        else
        if(opname == "XDn")
            {
            //m = -1 state along x axis
            Op = mixedIQTensor(dag(s),sP);
            Op = IQTensor(s);
            Op.set(Up,0.5);
            Op.set(Z0,-ISqrt2);
            Op.set(Dn,0.5);
            }
        else
        if(opname == "S2")
            {
            auto ssp1 = (s.m()==2 ? 0.75 : 2.);
            Op.set(Up,UpP,ssp1); 
            Op.set(Dn,DnP,ssp1);
            if(s.m() > 2)
                Op.set(Z0,Z0P,ssp1);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

inline SpinOne::
SpinOne(std::vector<IQIndex> const& inds)
    {
    int N = inds.size();
    auto sites = SiteStore(N);
    for(int j = 1, i = 0; j <= N; ++j, ++i)
        {
        auto& Ii = inds.at(i);
        if(Ii.m() != 3)
            {
            printfln("IQIndex at entry %d = %s",i,Ii);
            Error("Only S=1 IQIndices allowed in SpinOne(vector<IQIndex>) constructor");
            }
        sites.set(j,SpinOneSite(Ii));
        }
    SiteSet::init(std::move(sites));
    }

inline SpinOne::
SpinOne(int N, 
        Args const& args)
    {
    auto shedge = args.getBool("SHalfEdge",false);
    auto Lshedge = args.getBool("SHalfLeftEdge",false);

    auto sites = SiteStore(N);

    auto start = 1;
    if(shedge || Lshedge)
        {
        if(args.getBool("Verbose",false)) println("Placing a S=1/2 at site 1");
        sites.set(1,SpinHalfSite(1));
        start = 2;
        }

    for(int j = start; j < N; ++j)
        {
        sites.set(j,SpinOneSite(j));
        }

    if(shedge)
        {
        if(args.getBool("Verbose",false)) println("Placing a S=1/2 at site N=",N);
        sites.set(N,SpinHalfSite(N));
        }
    else
        {
        sites.set(N,SpinOneSite(N));
        }

    SiteSet::init(std::move(sites));
    }

void inline SpinOne::
read(std::istream& s)
    {
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        for(int j = 1; j <= N; ++j) 
            {
            auto I = IQIndex{};
            I.read(s);
            if(I.m() == 3) store.set(j,SpinOneSite(I));
            else if(I.m() == 2) store.set(j,SpinHalfSite(I));
            else Error(format("SpinOne cannot read index of size %d",I.m()));
            }
        init(std::move(store));
        }
    }

} //namespace itensor

#endif
