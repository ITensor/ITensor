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
        CHECK(1 == dim(i1));

        Index i2(1);
        CHECK(i2);
        CHECK(1 == dim(i2));
        }

    SECTION("Primes")
        {
        Index I(1);

        I = prime(I);
        CHECK_EQUAL(primeLevel(I),1);

        I = prime(I);
        CHECK_EQUAL(primeLevel(I),2);

        I = prime(I,7);
        CHECK_EQUAL(primeLevel(I),9);

        I = prime(I,-7);
        CHECK_EQUAL(primeLevel(I),2);

        I = setPrime(I,5);
        CHECK_EQUAL(primeLevel(I),5);
        
        I = noPrime(I);
        CHECK_EQUAL(primeLevel(I),0);

        }

    SECTION("IndexVal Basic")
        {
        Index i(4);
        IndexVal iv = i(2);
        CHECK(iv.val == 2);
        CHECK(iv.index == i);
        CHECK(dim(iv) == 4);
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
        CHECK(dim(i1) == dim(i));
        CHECK(tags(i1) == tags(i));
        CHECK(primeLevel(i1) == 0);

        auto i2 = sim(prime(i,3));
        CHECK(dim(i2) == dim(i));
        CHECK(tags(i2) == tags(i));
        CHECK(primeLevel(i2) == 0);
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

    SECTION("Tag Function")
        {
        auto i0 = Index(2);
        
        auto i = tags(i0, " -> i");
        CHECK(hasTags(i,"i"));

        i = tags(i, " -> site");
        CHECK(hasTags(i,"i,site"));

        i = tags(i, "i -> j");
        CHECK(hasTags(i,"j,site"));

        i = tags(i, "site -> ");
        CHECK(hasTags(i,"j"));

        i = tags(i, "j -> j,temp");
        CHECK(hasTags(i,"j,temp"));

        i = tags(i, "j,temp -> ");
        CHECK(i==i0);

        }
    }
