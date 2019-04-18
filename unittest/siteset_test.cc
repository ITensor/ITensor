#include "test.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/mps/sites/spintwo.h"
#include "itensor/mps/sites/electron.h"
#include "itensor/mps/sites/fermion.h"
#include "itensor/mps/sites/tj.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"

using namespace itensor;

TEST_CASE("SiteSetTest")
{

const int N = 10;

SECTION("Generic SiteSet")
    {
    auto sites = SiteSet(N,3);
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 3);
        CHECK(hasTags(sites(i),"Site,n="+str(i)));
        }
    }

SECTION("SpinHalf (QNs)")
    {
    auto sites = SpinHalf(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 2);
        CHECK(hasTags(sites(i),"Site,S=1/2,n="+str(i)));
        }

    op(sites,"Sz",2); 
    op(sites,"S+",2); 
    op(sites,"S-",2); 
    op(sites,"Sp",2); 
    op(sites,"Sm",2); 
    //TODO: test these throw the correct error
    //op(sites,"Sx",2); 
    //op(sites,"Sy",2); 
    //op(sites,"ISy",2); 
    }

SECTION("hasQNs")
    {
    auto sites = SpinHalf(N,{"ConserveQNs=",true});
    CHECK(hasQNs(sites));
    auto sitesNoQNs = SpinHalf(N,{"ConserveQNs=",false});
    CHECK(not hasQNs(sitesNoQNs));
    }

SECTION("SpinHalf (no QNs)")
    {
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 2);
        CHECK(hasTags(sites(i),"Site,S=1/2,n="+str(i)));
        }

    op(sites,"Sz",2); 
    op(sites,"S+",2); 
    op(sites,"S-",2); 
    op(sites,"Sp",2); 
    op(sites,"Sm",2); 
    op(sites,"Sx",2); 
    op(sites,"Sy",2); 
    op(sites,"ISy",2); 
    }

SECTION("SpinOne")
    {
    auto sites = SpinOne(N,{"ConserveQNs=",false});
    for(auto i : range1(N)) 
        {
        CHECK(dim(sites(i)) == 3);
        CHECK(hasTags(sites(i),"Site,S=1,n="+str(i)));
        }

    op(sites,"Sz",2); 
    op(sites,"S+",2); 
    op(sites,"S-",2); 
    op(sites,"Sp",2); 
    op(sites,"Sm",2); 
    op(sites,"Sx",2); 
    op(sites,"Sy",2); 
    op(sites,"ISy",2); 
    }

SECTION("SpinTwo")
    {
    auto sites = SpinTwo(N,{"ConserveQNs=",false});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 5);
        CHECK(hasTags(sites(i),"Site,S=2,n="+str(i)));
        }

    op(sites,"Sz",2);
    op(sites,"S+",2);
    op(sites,"S-",2);
    op(sites,"Sp",2);
    op(sites,"Sm",2);
    op(sites,"Sx",2);
    op(sites,"Sy",2);
    op(sites,"ISy",2);
    }

SECTION("Electron")
    {
    auto sites = Electron(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 4);
        CHECK(hasTags(sites(i),"Site"));
        }

    op(sites,"Nup",2); 
    op(sites,"Ndn",2); 
    op(sites,"Nupdn",2); 
    op(sites,"Ntot",2); 
    op(sites,"Sz",2); 
    op(sites,"Cup",2); 
    op(sites,"Cdn",2); 
    op(sites,"Aup",2); 
    op(sites,"Adn",2); 
    op(sites,"F",2); 
    }

SECTION("Fermion")
    {
    auto sites = Fermion(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 2);
        CHECK(hasTags(sites(i),"Site,Fermion,n="+str(i)));
        }

    op(sites,"N",2); 
    op(sites,"A",2); 
    op(sites,"Adag",2); 
    op(sites,"F",2); 
    }

SECTION("tJ")
    {
    auto sites = tJ(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 3);
        CHECK(hasTags(sites(i),"Site,tJ,n="+str(i)));
        }

    op(sites,"Nup",2); 
    op(sites,"Ndn",2); 
    op(sites,"Aup",2); 
    op(sites,"Adn",2); 
    op(sites,"F",2); 
    }
}

