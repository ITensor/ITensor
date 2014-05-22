#include "test.h"
#include "iqindex.h"

using namespace std;

TEST_CASE("IQIndexTest")
    {

SECTION("Null")
    {
    IQIndex i1;
    CHECK(i1.isNull());
    CHECK_EQUAL(1,i1.m());

    IQIndex I("I",Index("i"),QN());
    CHECK(!I.isNull());
    }

SECTION("Arrows")
    {
    CHECK_EQUAL(-In,Out);
    CHECK_EQUAL(-Out,In);
    }

SECTION("Primes")
    {
    IQIndex I("I",Index("i"),QN());

    I = primed(I);
    CHECK_EQUAL(I.primeLevel(),1);

    I = primed(I);
    CHECK_EQUAL(I.primeLevel(),2);

    I = primed(I,7);
    CHECK_EQUAL(I.primeLevel(),9);

    I = primed(I,-7);
    CHECK_EQUAL(I.primeLevel(),2);
    }


}
