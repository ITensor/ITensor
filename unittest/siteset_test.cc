#include "test.h"
#include "sites/spinhalf.h"
#include "sites/spinone.h"
#include "sites/hubbard.h"
#include "sites/spinless.h"
#include "sites/tj.h"

using namespace itensor;

TEST_CASE("SiteSetTest")
{

const int N = 10;

SECTION("SpinHalf")
    {
    SpinHalf sites(N);

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
    SpinOne sites(N);

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
    Hubbard sites(N);

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
    Spinless sites(N);

    sites.op("N",2); 
    sites.op("A",2); 
    sites.op("Adag",2); 
    sites.op("F",2); 
    }

SECTION("tJ")
    {
    tJ sites(N);

    sites.op("Nup",2); 
    sites.op("Ndn",2); 
    sites.op("Aup",2); 
    sites.op("Adn",2); 
    sites.op("F",2); 
    }
}

