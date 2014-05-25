#include "test.h"
#include "index.h"

using namespace itensor;
using namespace std;


TEST_CASE("IndexTest")
    {
    SECTION("Boolean")
        {
        Index i1;
        CHECK(!i1);
        CHECK(1 == i1.m());

        Index i2("i2");
        CHECK(i2);
        CHECK(1 == i2.m());
        }

    SECTION("Primes")
        {
        Index I("I");

        I = prime(I);
        CHECK_EQUAL(I.primeLevel(),1);

        I = prime(I);
        CHECK_EQUAL(I.primeLevel(),2);

        I = prime(I,7);
        CHECK_EQUAL(I.primeLevel(),9);

        I = prime(I,-7);
        CHECK_EQUAL(I.primeLevel(),2);
        }
    }


