#include "test.h"
#include "index.h"
#include <boost/test/unit_test.hpp>

using namespace itensor;
using namespace std;

BOOST_AUTO_TEST_SUITE(IndexTest)

TEST(Null)
    {
    Index i1;
    CHECK(i1.isNull());
    CHECK_EQUAL(1,i1.m());

    Index i2("i2");
    CHECK(!i2.isNull());
    }

TEST(Primes)
    {
    Index I("I");

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

