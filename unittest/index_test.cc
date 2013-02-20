#include "test.h"
#include "index.h"
#include <boost/test/unit_test.hpp>

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

TEST(Complex)
    {
    Index I = Index::IndReIm();

    I.prime(All);
    CHECK_EQUAL(I.primeLevel(),0);

    I.prime();
    CHECK_EQUAL(I.primeLevel(),1);

    I.prime(2);
    CHECK_EQUAL(I.primeLevel(),3);

    Index J = primed(Index::IndReIm());
    CHECK_EQUAL(J.primeLevel(),1);

    Index K = primed(Index::IndReIm(),All);
    CHECK_EQUAL(K.primeLevel(),0);

    Index L = Index::IndReIm();
    L.prime(All);
    CHECK_EQUAL(L.primeLevel(),0);

    CHECK_EQUAL(Index::IndReIm().primeLevel(),0);
    CHECK_EQUAL(Index::IndReImP().primeLevel(),1);
    CHECK_EQUAL(Index::IndReImPP().primeLevel(),2);
    }

BOOST_AUTO_TEST_SUITE_END()

