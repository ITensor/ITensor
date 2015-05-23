#include "test.h"
#include "itensor/iqindex.h"

using namespace itensor;
using namespace std;

TEST_CASE("IQIndexTest")
    {

SECTION("Null")
    {
    IQIndex i1;
    CHECK(!i1);
    CHECK_EQUAL(1,i1.m());

    IQIndex I("I",Index("i"),QN());
    CHECK(I);
    }

SECTION("Arrows")
    {
    CHECK_EQUAL(-In,Out);
    CHECK_EQUAL(-Out,In);
    }

SECTION("Primes")
    {
    IQIndex I("I",Index("i"),QN());

    I = prime(I);
    CHECK_EQUAL(I.primeLevel(),1);

    I = prime(I);
    CHECK_EQUAL(I.primeLevel(),2);

    I = prime(I,7);
    CHECK_EQUAL(I.primeLevel(),9);

    I = prime(I,-7);
    CHECK_EQUAL(I.primeLevel(),2);
    }

SECTION("Constructors")
    {
    SECTION("One Index")
        {
        auto i1 = Index("i1",5,Site);
        auto I = IQIndex("one",
                         i1,QN(-1));
        CHECK(I.nindex() == 1);
        CHECK(I[0] == i1);
        }

    SECTION("Two Indices")
        {
        auto i1 = Index("i1",5,Site);
        auto i2 = Index("i2",8,Site);
        auto I = IQIndex("two",
                         i1,QN(-1),
                         i2,QN(0));
        CHECK(I.nindex() == 2);
        CHECK(I[0] == i1);
        CHECK(I[1] == i2);
        }

    SECTION("Three Indices")
        {
        auto i1 = Index("i1",5,Site);
        auto i2 = Index("i2",8,Site);
        auto i3 = Index("i3",2,Site);
        auto I = IQIndex("three",
                         i1,QN(-1),
                         i2,QN(0),
                         i3,QN(+1));
        CHECK(I.nindex() == 3);
        CHECK(I[0] == i1);
        CHECK(I[1] == i2);
        CHECK(I[2] == i3);
        CHECK(I.qn(1) == QN(-1));
        CHECK(I.qn(2) == QN(0));
        CHECK(I.qn(3) == QN(+1));
        CHECK(I.dir() == Out);
        }

    SECTION("Four Indices")
        {
        auto i1 = Index("i1",5);
        auto i2 = Index("i2",8);
        auto i3 = Index("i3",2);
        auto i4 = Index("i4",4);
        auto I = IQIndex("four",
                         In,
                         i1,QN(-1),
                         i2,QN(0),
                         i3,QN(+1),
                         i4,QN(+2));
        CHECK(I.nindex() == 4);
        CHECK(I.dir() == In);
        CHECK(I[0] == i1);
        CHECK(I[1] == i2);
        CHECK(I[2] == i3);
        CHECK(I[3] == i4);
        CHECK(I.qn(1) == QN(-1));
        CHECK(I.qn(2) == QN(0));
        CHECK(I.qn(3) == QN(+1));
        CHECK(I.qn(4) == QN(+2));
        }
    }


}
