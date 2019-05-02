#pragma once

#include "itensor/mps/mpo.h"
#include "itensor/util/print_macro.h"

namespace itensor {

MPO
makeS2(SiteSet const& sites,
       Args const& args = Args::global())
    {
    auto skip_ancilla = args.getBool("SkipAncilla",false);
    auto N = length(sites);

    auto S2 = MPO(sites);

    auto links = std::vector<Index>(N+1);
    for(auto n : range(N+1))
        {
        auto ts = format("Link,l=%d",n);
        links.at(n) = Index(QN({"Sz",0}),3,
                            QN({"Sz",-2}),1,
                            QN({"Sz",+2}),1,ts);
        }

    for(auto n : range1(N))
        {
        bool phys_site = true;
        if(skip_ancilla && (n%2==0)) phys_site = false;

        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = S2.ref(n);
        W = ITensor(row,col,dag(sites(n)),prime(sites(n)));

        W += op(sites,"Id",n) * setElt(row(1)) * setElt(col(1));
        W += op(sites,"Id",n) * setElt(row(2)) * setElt(col(2));

        if(phys_site)
            {
            W += op(sites,"S2",n) * setElt(row(2)) * setElt(col(1));
            }

        if(phys_site)
            {
            W += 2*op(sites,"Sz",n) * setElt(row(2)) * setElt(col(3));
            W +=   op(sites,"Sz",n) * setElt(row(3)) * setElt(col(1));
            }
        W += op(sites,"Id",n) * setElt(row(3)) * setElt(col(3));

        if(phys_site)
            {
            W += op(sites,"S+",n) * setElt(row(2)) * setElt(col(4));
            W += op(sites,"S-",n) * setElt(row(4)) * setElt(col(1));
            }
        W += op(sites,"Id",n) * setElt(row(4)) * setElt(col(4));

        if(phys_site)
            {
            W += op(sites,"S-",n) * setElt(row(2)) * setElt(col(5));
            W += op(sites,"S+",n) * setElt(row(5)) * setElt(col(1));
            }
        W += op(sites,"Id",n) * setElt(row(5)) * setElt(col(5));

        //W.scaleTo(1.);
        }

    S2.ref(1) *= setElt(links.at(0)(2));
    S2.ref(N) *= setElt(dag(links.at(N))(1));

    return S2;
    }

MPO
makeTotSz2(SiteSet const& sites,
           Args const& args = Args::global())
    {
    auto skip_ancilla = args.getBool("SkipAncilla",false);
    auto N = sites.N();

    auto Sz2 = MPO(sites);

    auto links = std::vector<Index>(N+1);
    for(auto n : range(N+1))
        {
        auto ts = format("Link,MPO,%d",n);
        links.at(n) = Index(QN({"Sz=",0}),3,ts);
        }

    for(auto n : range1(N))
        {
        bool phys_site = true;
        if(skip_ancilla && (n%2==0)) phys_site = false;

        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = Sz2.ref(n);
        W = ITensor(row,col,dag(sites(n)),prime(sites(n)));

        W += op(sites,"Id",n) * setElt(row(1)) * setElt(col(1));
        W += op(sites,"Id",n) * setElt(row(2)) * setElt(col(2));

        if(phys_site)
            {
            W += op(sites,"Sz*Sz",n) * setElt(row(2)) * setElt(col(1));
            }

        if(phys_site)
            {
            W += 2*op(sites,"Sz",n) * setElt(row(2)) * setElt(col(3));
            W +=   op(sites,"Sz",n) * setElt(row(3)) * setElt(col(1));
            }
        W += op(sites,"Id",n) * setElt(row(3)) * setElt(col(3));

        //W.scaleTo(1.);
        }

    Sz2.ref(1) *= setElt(links.at(0)(2));
    Sz2.ref(N) *= setElt(dag(links.at(N))(1));

    return Sz2;
    }

MPO
makeSxy2(SiteSet const& sites,
         Args const& args = Args::global())
    {
    auto skip_ancilla = args.getBool("SkipAncilla",false);
    auto N = sites.N();

    auto Sxy2 = MPO(sites);

    auto links = std::vector<Index>(N+1);
    for(auto n : range(N+1))
        {
        auto ts = format("Link,MPO,%d",n);
        links.at(n) = Index(QN({"Sz",0}),3,
                            QN({"Sz",-2}),1,
                            QN({"Sz",+2}),1,ts);
        }

    for(auto n : range1(N))
        {
        bool phys_site = true;
        if(skip_ancilla && (n%2==0)) phys_site = false;

        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = Sxy2.ref(n);
        W = ITensor(row,col,dag(sites(n)),prime(sites(n)));

        W += op(sites,"Id",n) * setElt(row(1)) * setElt(col(1));
        W += op(sites,"Id",n) * setElt(row(2)) * setElt(col(2));

        if(phys_site)
            {
            W += 0.5*op(sites,"Id",n) * setElt(row(2)) * setElt(col(1));
            }

        //if(phys_site)
        //    {
        //    W += 2*op(sites,"Sz",n) * setElt(row(2)) * setElt(col(3));
        //    W +=   op(sites,"Sz",n) * setElt(row(3)) * setElt(col(1));
        //    }
        //W += op(sites,"Id",n) * setElt(row(3)) * setElt(col(3));

        if(phys_site)
            {
            W += op(sites,"S+",n) * setElt(row(2)) * setElt(col(4));
            W += op(sites,"S-",n) * setElt(row(4)) * setElt(col(1));
            }
        W += op(sites,"Id",n) * setElt(row(4)) * setElt(col(4));

        if(phys_site)
            {
            W += op(sites,"S-",n) * setElt(row(2)) * setElt(col(5));
            W += op(sites,"S+",n) * setElt(row(5)) * setElt(col(1));
            }
        W += op(sites,"Id",n) * setElt(row(5)) * setElt(col(5));

        //W.scaleTo(1.);
        }

    Sxy2.ref(1) *= setElt(links.at(0)(2));
    Sxy2.ref(N) *= setElt(dag(links.at(N))(1));

    return Sxy2;
    }

} //namespace itensor
