#pragma once

#include "itensor/mps/mpo.h"
#include "itensor/util/print_macro.h"

namespace itensor {

IQMPO
makeS2(SiteSet const& sites,
       Args const& args = Args::global())
    {
    auto skip_ancilla = args.getBool("SkipAncilla",false);
    auto N = sites.N();

    auto S2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Sz=",0),
                              Index("+",1),QN("Sz=",-2),
                              Index("-",1),QN("Sz=",+2));
        }

    for(auto n : range1(N))
        {
        bool phys_site = true;
        if(skip_ancilla && (n%2==0)) phys_site = false;

        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = S2.Aref(n);
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));

        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);

        if(phys_site)
            {
            W += sites.op("S2",n) * row(2) * col(1);
            }

        if(phys_site)
            {
            W += 2*sites.op("Sz",n) * row(2) * col(3);
            W +=   sites.op("Sz",n) * row(3) * col(1);
            }
        W += sites.op("Id",n) * row(3) * col(3);

        if(phys_site)
            {
            W += sites.op("S+",n) * row(2) * col(4);
            W += sites.op("S-",n) * row(4) * col(1);
            }
        W += sites.op("Id",n) * row(4) * col(4);

        if(phys_site)
            {
            W += sites.op("S-",n) * row(2) * col(5);
            W += sites.op("S+",n) * row(5) * col(1);
            }
        W += sites.op("Id",n) * row(5) * col(5);

        //W.scaleTo(1.);
        }

    S2.Aref(1) *= setElt(links.at(0)(2));
    S2.Aref(N) *= setElt(dag(links.at(N))(1));

    return S2;
    }

IQMPO
makeTotSz2(SiteSet const& sites,
           Args const& args = Args::global())
    {
    auto skip_ancilla = args.getBool("SkipAncilla",false);
    auto N = sites.N();

    auto Sz2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Sz=",0));
        }

    for(auto n : range1(N))
        {
        bool phys_site = true;
        if(skip_ancilla && (n%2==0)) phys_site = false;

        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = Sz2.Aref(n);
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));

        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);

        if(phys_site)
            {
            W += sites.op("Sz*Sz",n) * row(2) * col(1);
            }

        if(phys_site)
            {
            W += 2*sites.op("Sz",n) * row(2) * col(3);
            W +=   sites.op("Sz",n) * row(3) * col(1);
            }
        W += sites.op("Id",n) * row(3) * col(3);

        //W.scaleTo(1.);
        }

    Sz2.Aref(1) *= setElt(links.at(0)(2));
    Sz2.Aref(N) *= setElt(dag(links.at(N))(1));

    return Sz2;
    }

IQMPO
makeSxy2(SiteSet const& sites,
         Args const& args = Args::global())
    {
    auto skip_ancilla = args.getBool("SkipAncilla",false);
    auto N = sites.N();

    auto Sxy2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Sz=",0),
                              Index("+",1),QN("Sz=",-2),
                              Index("-",1),QN("Sz=",+2));
        }

    for(auto n : range1(N))
        {
        bool phys_site = true;
        if(skip_ancilla && (n%2==0)) phys_site = false;

        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = Sxy2.Aref(n);
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));

        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);

        if(phys_site)
            {
            W += 0.5*sites.op("Id",n) * row(2) * col(1);
            }

        //if(phys_site)
        //    {
        //    W += 2*sites.op("Sz",n) * row(2) * col(3);
        //    W +=   sites.op("Sz",n) * row(3) * col(1);
        //    }
        //W += sites.op("Id",n) * row(3) * col(3);

        if(phys_site)
            {
            W += sites.op("S+",n) * row(2) * col(4);
            W += sites.op("S-",n) * row(4) * col(1);
            }
        W += sites.op("Id",n) * row(4) * col(4);

        if(phys_site)
            {
            W += sites.op("S-",n) * row(2) * col(5);
            W += sites.op("S+",n) * row(5) * col(1);
            }
        W += sites.op("Id",n) * row(5) * col(5);

        //W.scaleTo(1.);
        }

    Sxy2.Aref(1) *= setElt(links.at(0)(2));
    Sxy2.Aref(N) *= setElt(dag(links.at(N))(1));

    return Sxy2;
    }

} //namespace itensor
