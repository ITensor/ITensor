//
// Distributed under the ITensor Library License, Version 1.2
// (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINTWO_H
#define __ITENSOR_SPINTWO_H
#include "itensor/mps/siteset.h"
#include "itensor/mps/sites/spinhalf.h"

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
    IQIndex s;
	public:

    SpinTwoSite() { }

    SpinTwoSite(IQIndex I) : s(I) { }

    SpinTwoSite(int n, Args const& args = Args::global())
		{
        s = IQIndex{nameint("S=2 site=",n),
            Index(nameint("Up:site",n),1,Site),QN("Sz=",+4),
            Index(nameint("Upi:site",n),1,Site),QN("Sz=",+2),
            Index(nameint("Z0:site",n),1,Site),QN("Sz=",0),
            Index(nameint("Dni:site",n),1,Site),QN("Sz=",-2),
            Index(nameint("Dn:site",n),1,Site),QN("Sz=",-4)};
		}

    IQIndex
    index() const { return s; }

    IQIndexVal
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
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
		}

    IQTensor
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

        auto Op = IQTensor(dag(s),sP);

        if (opname == "Sz")
            {
            Op.set(Up,UpP,+2.0);
            Op.set(Upi,UpiP,+1.0);
            Op.set(Dni,DniP,-1.0);
            Op.set(Dn,DnP,-2.0);
            }
        else if (opname == "Sx")
            {
            //mixedIQTensor call needed here
            //because as an IQTensor, Op would
            //not have a well defined QN flux
            Op = mixedIQTensor(s,sP);
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
            //mixedIQTensor call needed here
            //because as an IQTensor, Op would
            //not have a well defined QN flux
            Op = mixedIQTensor(s,sP);
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
            //mixedIQTensor call needed here
            //because as an IQTensor, Op would
            //not have a well defined QN flux
            Op = mixedIQTensor(s,sP);
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
            Error("Operator " + opname + " name not recognized");
            }

        return Op;
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
        sites.set(1,SpinHalfSite(1));
        start = 2;
        }

    for(int j = start; j < N; ++j)
        {
        sites.set(j,SpinTwoSite(j));
        }

    if(shedge)
        {
        if(args.getBool("Verbose",false)) println("Placing a S=1/2 at site N=",N);
        sites.set(N,SpinHalfSite(N));
        }
    else
        {
        sites.set(N,SpinTwoSite(N));
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
            auto I = IQIndex{};
            I.read(s);
            if(I.m() == 3) store.set(j,SpinTwoSite(I));
            else if(I.m() == 2) store.set(j,SpinHalfSite(I));
            else Error(format("SpinTwo cannot read index of size %d",I.m()));
            }
        init(std::move(store));
        }
	}

} //namespace itensor

#endif
