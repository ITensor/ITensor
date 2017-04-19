#include "test.h"
#include "itensor/index.h"
#include "itensor/util/print_macro.h"

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

    SECTION("IndexVal Basic")
        {
        Index i("i",4);
        IndexVal iv = i(2);
        CHECK(iv.val == 2);
        CHECK(iv.index == i);
        CHECK(iv.m() == 4);
        if(iv) CHECK(true);
        else   CHECK(false);

        IndexVal ivP = prime(iv);
        iv.prime();
        CHECK(ivP == iv);

        IndexVal def;
        if(def) CHECK(false);
        else    CHECK(true);

        }

    SECTION("sim function")
        {
        auto i = Index("i",4,Atype);

        auto i1 = sim(i);
        CHECK(i1.m() == i.m());
        CHECK(i1.type() == i.type());
        CHECK(i1.primeLevel() == 0);

        auto i2 = sim(prime(i,3));
        CHECK(i2.m() == i.m());
        CHECK(i2.type() == i.type());
        CHECK(i2.primeLevel() == 0);
        }
    }
