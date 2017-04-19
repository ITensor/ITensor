#include "test.h"
#include "itensor/iqindex.h"
#include "itensor/util/print_macro.h"

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
                         i1,QN(-1),
                         i2,QN(0),
                         i3,QN(+1),
                         i4,QN(+2),
                         In);
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

SECTION("Iterator")
    {
    auto is = vector<Index>(4+1);
    is[1] = Index("i1",1);
    is[2] = Index("i2",2);
    is[3] = Index("i3",3);
    is[4] = Index("i4",4);

    auto I = IQIndex("I",
                     is[1],QN(-1),
                     is[2],QN(0),
                     is[3],QN(+1),
                     is[4],QN(+2));

    auto n = 1;
    for(auto i : I)
        {
        CHECK(i.index == is[n]);
        ++n;
        }
    CHECK(n == 1+I.nindex());

    auto I4 = prime(I,4);
    n = 1;
    for(auto i : I4)
        {
        CHECK(i.index == prime(is[n],4));
        ++n;
        }
    CHECK(n == 1+I.nindex());
    }

SECTION("sim function")
    {
    auto i1 = Index("i1",5,Site);
    auto i2 = Index("i2",8,Site);
    auto i3 = Index("i3",2,Site);
    auto I = IQIndex("I",
                     i1,QN(-1),
                     i2,QN(0),
                     i3,QN(+1));

    auto S1 = sim(I);
    for(auto n : range1(I.nindex()))
        {
        CHECK(I.index(n) == S1.index(n));
        }
    CHECK(S1 != I);
    CHECK(S1.type() == I.type());
    CHECK(S1.m() == I.m());
    CHECK(S1.dir() == I.dir());
    CHECK(S1.primeLevel() == 0);

    auto S2 = sim(prime(I));
    for(auto n : range1(I.nindex()))
        {
        CHECK(I.index(n) == S2.index(n));
        }
    CHECK(S2 != I);
    CHECK(S2.type() == I.type());
    CHECK(S2.m() == I.m());
    CHECK(S2.dir() == I.dir());
    CHECK(S2.primeLevel() == 0);
    }


}
