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

TEST(Complex)
    {
    IQIndex I = IQIndex::IndReIm();

    I.prime(All);
    CHECK_EQUAL(I.primeLevel(),0);

    I.prime();
    CHECK_EQUAL(I.primeLevel(),1);

    I.prime(2);
    CHECK_EQUAL(I.primeLevel(),3);

    IQIndex J = primed(IQIndex::IndReIm());
    CHECK_EQUAL(J.primeLevel(),1);

    IQIndex K = primed(IQIndex::IndReIm(),All);
    CHECK_EQUAL(K.primeLevel(),0);

    IQIndex L = IQIndex::IndReIm();
    L.prime(All);
    CHECK_EQUAL(L.primeLevel(),0);

    CHECK_EQUAL(IQIndex::IndReIm().primeLevel(),0);
    CHECK_EQUAL(IQIndex::IndReImP().primeLevel(),1);
    CHECK_EQUAL(IQIndex::IndReImPP().primeLevel(),2);
    }

BOOST_AUTO_TEST_SUITE_END()

