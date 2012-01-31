#include "test.h"
#include "index.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(IndexTest)

BOOST_AUTO_TEST_CASE(Null)
{
    Index i1;
    CHECK(i1.isNull());
    CHECK_EQUAL(1,i1.m());

    Index i2("i2");
    CHECK(i2.isNotNull());
}

BOOST_AUTO_TEST_CASE(Primes)
{
    Index I("I");

    I = I.primed();
    CHECK_EQUAL(I.primeLevel(),1);

    I = I.primed();
    CHECK_EQUAL(I.primeLevel(),2);

    I = I.primed(7);
    CHECK_EQUAL(I.primeLevel(),9);

    I = I.primed(-7);
    CHECK_EQUAL(I.primeLevel(),2);

}

BOOST_AUTO_TEST_SUITE_END()

