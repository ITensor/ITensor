#include "test.h"
#include "combiner.h"
#include <boost/test/unit_test.hpp>

using namespace itensor;

struct CombinerDefaults
{
    const Index s1,s2,s3,s4,
          s1P,s2P,s3P,s4P,
          l1,l2,l3,l4,l5,l6,l7,l8,
          a1,a2,a3,a4,
          b2,b3,b4,b5;

    CombinerDefaults() :
    s1(Index("s1",2,Site)),
    s2(Index("s2",2,Site)),
    s3(Index("s3",2,Site)),
    s4(Index("s4",2,Site)),
    s1P(prime(s1)),
    s2P(prime(s2)),
    s3P(prime(s3)),
    s4P(prime(s4)),
    l1(Index("l1",2)),
    l2(Index("l2",2)),
    l3(Index("l3",2)),
    l4(Index("l4",2)),
    l5(Index("l5",2)),
    l6(Index("l6",2)),
    l7(Index("l7",2)),
    l8(Index("l8",2)),
    a1(Index("a1")),
    a2(Index("a2")),
    a3(Index("a3")),
    a4(Index("a4")),
    b2(Index("b2",2)),
    b3(Index("b3",3)),
    b4(Index("b4",4)),
    b5(Index("b5",5))
    {

    }

    ~CombinerDefaults() { }

};

BOOST_FIXTURE_TEST_SUITE(CombinerTest,CombinerDefaults)

TEST(Constructors)
{

    Combiner c1;

    Combiner c2(l1,a1,l2);

    c2.init();

    CHECK(c2.isInit());
    CHECK(hasindex(c2,l1));
    CHECK(hasindex(c2,a1));
    CHECK(hasindex(c2,l2));
    CHECK_EQUAL(c2.right().m(),l1.m()*a1.m()*l2.m());

}

TEST(addLeft)
{
    Combiner c1;

    c1.addleft(l2);
    c1.addleft(a2);
    c1.addleft(l3);
    c1.addleft(l4);

    CHECK(hasindex(c1,l2));
    CHECK(hasindex(c1,l3));
    CHECK(hasindex(c1,a2));
    CHECK(hasindex(c1,l4));

    c1.init("cname");

    CHECK_EQUAL(c1.right().name(),"cname");
    CHECK_EQUAL(c1.right().m(),l2.m()*l3.m()*l4.m()*a2.m());

}

TEST(Product)
{

    ITensor A(a1,b3,l2,a4,l3);
    A.randomize();

    Combiner c; 
    c.addleft(b3);
    c.addleft(l3);
    c.addleft(a1);

    c.init();
    Index r = c.right();

    ITensor cA = c * A;

    CHECK(hasindex(cA,l2));
    CHECK(hasindex(cA,a4));
    CHECK(hasindex(cA,r));
    CHECK(!hasindex(cA,b3));
    CHECK(!hasindex(cA,l3));
    CHECK(!hasindex(cA,a1));

    CHECK_CLOSE(A(l2(1),l3(1),b3(1)),cA(l2(1),r(1)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(1),b3(2)),cA(l2(1),r(2)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(1),b3(3)),cA(l2(1),r(3)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(2),b3(1)),cA(l2(1),r(4)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(2),b3(2)),cA(l2(1),r(5)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(2),b3(3)),cA(l2(1),r(6)),1E-10);

    ITensor ucA = c * cA;

    CHECK(hasindex(ucA,a1));
    CHECK(hasindex(ucA,b3));
    CHECK(hasindex(ucA,l2));
    CHECK(hasindex(ucA,a4));
    CHECK(hasindex(ucA,l3));
    CHECK(!hasindex(ucA,r));

    for(int k2 = 1; k2 <= 2; ++k2)
    for(int k3 = 1; k3 <= 2; ++k3)
    for(int j3 = 1; j3 <= 2; ++j3)
    {
        CHECK_CLOSE(ucA(l2(k2),l3(k3),b3(j3)),A(l2(k2),l3(k3),b3(j3)),1E-10);
    }

    ITensor B(a1,b3,a3,l2,a4,l3);
    B.randomize();

    Combiner c2; 
    c2.addleft(a3);
    c2.addleft(a1);

    ITensor cB = c2 * B;

    Index r2 = c2.right();

    CHECK(hasindex(cB,b3));
    CHECK(hasindex(cB,l2));
    CHECK(hasindex(cB,a4));
    CHECK(hasindex(cB,l3));
    CHECK(!hasindex(cB,a3));
    CHECK(!hasindex(cB,a1));

    for(int k2 = 1; k2 <= 2; ++k2)
    for(int k3 = 1; k3 <= 2; ++k3)
    for(int j3 = 1; j3 <= 2; ++j3)
    {
        CHECK_CLOSE(cB(b3(j3),l2(k2),l3(k3)),B(b3(j3),l2(k2),l3(k3)),1E-10);
    }

}

BOOST_AUTO_TEST_SUITE_END()
