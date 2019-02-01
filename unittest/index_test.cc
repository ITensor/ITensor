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

        Index i2(1);
        CHECK(i2);
        CHECK(1 == i2.m());
        }

    SECTION("Primes")
        {
        Index I(1);

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
        Index i(4);
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
        auto i = Index(4,"i");

        auto i1 = sim(i);
        CHECK(i1.m() == i.m());
        CHECK(i1.tags() == i.tags());
        CHECK(i1.primeLevel() == 0);

        auto i2 = sim(prime(i,3));
        CHECK(i2.m() == i.m());
        CHECK(i2.tags() == i.tags());
        CHECK(i2.primeLevel() == 0);
        }

    SECTION("Tag Basics")
        {
        auto i = Index(3);

        auto ia = addTags(i,"a");
        CHECK(hasTags(ia,"a"));
        CHECK(ia != i);
        CHECK(ia == ia);

        auto ib = addTags(i,"b");
        CHECK(hasTags(ib,"b"));
        CHECK(ib != i);
        CHECK(ib == ib);

        CHECK(ia != ib);

        auto iab1 = addTags(ia,"b");
        auto iab2 = addTags(ib,"a");
        CHECK(iab1 == iab2);
        CHECK(hasTags(iab1,"a"));
        CHECK(hasTags(iab1,"b"));
        CHECK(hasTags(iab1,"a,b"));

        CHECK(removeTags(iab1,"a,b") == i);
        CHECK(removeTags(iab1,"b") == ia);
        CHECK(removeTags(iab1,"a") == ib);

        auto ic = setTags(iab1,"c");
        CHECK(hasTags(ic,"c"));

        CHECK(hasTags(replaceTags(ic,"c","a"),"a"));
        CHECK(hasTags(replaceTags(ic,"c","a,b"),"a,b"));
        CHECK(hasTags(replaceTags(ic,"","a,b"),"a,b"));
        CHECK(!hasTags(replaceTags(ic,"c","a,b"),"c"));
        CHECK(hasTags(addTags(ic,"a,b"),"a,b,c"));
        }

    }
