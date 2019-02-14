#include "test.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/mps/sites/hubbard.h"
#include "itensor/mps/sites/spinless.h"
#include "itensor/mps/sites/tj.h"
#include "itensor/util/print_macro.h"

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
        CHECK(hasTags(sites(i),"Site"));
        }
    }

SECTION("SpinHalf (QNs)")
    {
    auto sites = SpinHalf(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 2);
        CHECK(hasTags(sites(i),"Site"));
        }

    sites.op("Sz",2); 
    sites.op("S+",2); 
    sites.op("S-",2); 
    sites.op("Sp",2); 
    sites.op("Sm",2); 
    //TODO: test these throw the correct error
    //sites.op("Sx",2); 
    //sites.op("Sy",2); 
    //sites.op("ISy",2); 
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
        CHECK(hasTags(sites(i),"Site"));
        }

    sites.op("Sz",2); 
    sites.op("S+",2); 
    sites.op("S-",2); 
    sites.op("Sp",2); 
    sites.op("Sm",2); 
    sites.op("Sx",2); 
    sites.op("Sy",2); 
    sites.op("ISy",2); 
    }

SECTION("SpinOne")
    {
    auto sites = SpinOne(N,{"ConserveQNs=",false});
    for(auto i : range1(N)) 
        {
        CHECK(dim(sites(i)) == 3);
        CHECK(hasTags(sites(i),"Site"));
        }

    sites.op("Sz",2); 
    sites.op("S+",2); 
    sites.op("S-",2); 
    sites.op("Sp",2); 
    sites.op("Sm",2); 
    sites.op("Sx",2); 
    sites.op("Sy",2); 
    sites.op("ISy",2); 
    }

SECTION("Hubbard")
    {
    auto sites = Hubbard(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 4);
        CHECK(hasTags(sites(i),"Site"));
        }

    sites.op("Nup",2); 
    sites.op("Ndn",2); 
    sites.op("Nupdn",2); 
    sites.op("Ntot",2); 
    sites.op("Sz",2); 
    sites.op("Cup",2); 
    sites.op("Cdn",2); 
    sites.op("Aup",2); 
    sites.op("Adn",2); 
    sites.op("F",2); 
    }

SECTION("Spinless")
    {
    auto sites = Spinless(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(sites(i).dim() == 2);
        CHECK(hasTags(sites(i),"Site"));
        }

    sites.op("N",2); 
    sites.op("A",2); 
    sites.op("Adag",2); 
    sites.op("F",2); 
    }

SECTION("tJ")
    {
    auto sites = tJ(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(sites(i).dim() == 3);
        CHECK(hasTags(sites(i),"Site"));
        }

    sites.op("Nup",2); 
    sites.op("Ndn",2); 
    sites.op("Aup",2); 
    sites.op("Adn",2); 
    sites.op("F",2); 
    }
}

