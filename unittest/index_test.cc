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
        CHECK(hasTags(i2,"i"));
        CHECK(primeLevel(i2) == 3);
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

        auto ic = setTags(iab1,"c,0");
        CHECK(hasTags(ic,"c"));

        CHECK(ic == setTags(iab1,"c"));
        CHECK(i == noTags(ic));
        CHECK(noTags(ia) == noTags(ic));
        CHECK(hasTags(replaceTags(ic,"c","a"),"a"));
        CHECK(hasTags(replaceTags(ic,"c","a,b"),"a,b"));
        CHECK(hasTags(replaceTags(ic,"","a,b"),"a,b"));
        CHECK(!hasTags(replaceTags(ic,"c","a,b"),"c"));
        CHECK(hasTags(addTags(ic,"a,b"),"a,b,c"));
        }

    SECTION("Ignore spaces input string tests")
      {
      auto ts1 = TagSet("a,n=1");
      auto ts2 = TagSet(" a, n = 1 ");
      CHECK(ts1 == ts2);

      ts1.addTags(" b ");
      CHECK(ts1 == TagSet("a, b, n = 1"));
      }

    SECTION("Access Prime Level from TagSet")
      {
      auto s = Index(3);

      auto i = addTags(s,"i");
      CHECK(hasTags(i,"i"));
      CHECK(hasTags(i,"0"));
      CHECK(primeLevel(i)==0);

      auto ip = replaceTags(i,"0","1");
      CHECK(hasTags(ip,"i"));
      CHECK(hasTags(ip,"1"));
      CHECK(primeLevel(ip)==1);
      CHECK(ip == replaceTags(s,"0","i,1"));

      auto x = replaceTags(ip,"i,1","x,3");
      CHECK(hasTags(x,"x"));
      CHECK(hasTags(x,"x,3"));
      CHECK(primeLevel(x)==3);
      CHECK(x == replaceTags(s,"0","x,3"));


      }

    }
