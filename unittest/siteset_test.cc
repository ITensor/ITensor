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
        CHECK(sites(i).m() == 3);
        CHECK(sites(i).type() == Site);
        }
    }

SECTION("SpinHalf")
    {
    auto sites = SpinHalf(N);
    for(auto i : range1(N))
        {
        CHECK(sites(i).m() == 2);
        CHECK(sites(i).type() == Site);
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
    auto sites = SpinOne(N);
    for(auto i : range1(N)) 
        {
        CHECK(sites(i).m() == 3);
        CHECK(sites(i).type() == Site);
        }

    sites.op("Sz",2); 
    sites.op("S+",2); 
    sites.op("S-",2); 
    sites.op("Sp",2); 
    sites.op("Sm",2); 
    sites.op("Sx",2); 
    //sites.op("Sy",2); 
    sites.op("ISy",2); 
    }

SECTION("Hubbard")
    {
    auto sites = Hubbard(N);
    for(auto i : range1(N))
        {
        CHECK(sites(i).m() == 4);
        CHECK(sites(i).type() == Site);
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
    auto sites = Spinless(N);
    for(auto i : range1(N))
        {
        CHECK(sites(i).m() == 2);
        CHECK(sites(i).type() == Site);
        }

    sites.op("N",2); 
    sites.op("A",2); 
    sites.op("Adag",2); 
    sites.op("F",2); 
    }

SECTION("tJ")
    {
    auto sites = tJ(N);
    for(auto i : range1(N))
        {
        CHECK(sites(i).m() == 3);
        CHECK(sites(i).type() == Site);
        }

    sites.op("Nup",2); 
    sites.op("Ndn",2); 
    sites.op("Aup",2); 
    sites.op("Adn",2); 
    sites.op("F",2); 
    }
}

