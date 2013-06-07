#include "test.h"
#include "iqindex.h"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_SUITE(IQIndexTest)

TEST(Null)
    {
    IQIndex i1;
    CHECK(i1.isNull());
    CHECK_EQUAL(1,i1.m());

    IQIndex I("I",Index("i"),QN());
    CHECK(!I.isNull());
    }

TEST(Arrows)
    {
    CHECK_EQUAL(-In,Out);
    CHECK_EQUAL(-Out,In);
    }

TEST(Primes)
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


BOOST_AUTO_TEST_SUITE_END()

