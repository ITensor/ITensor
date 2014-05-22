#include "test.h"
#include "index.h"

using namespace std;


TEST_CASE("Null")
    {
    Index i1;
    CHECK(i1.isNull());
    CHECK_EQUAL(1,i1.m());

    Index i2("i2");
    CHECK(!i2.isNull());
    }

TEST_CASE("Primes")
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



