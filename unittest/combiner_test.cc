#include "test.h"
#include "combiner.h"
#include <boost/test/unit_test.hpp>

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
    s1P(s1.primed()),
    s2P(s2.primed()),
    s3P(s3.primed()),
    s4P(s4.primed()),
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

BOOST_AUTO_TEST_CASE(Constructors)
{

    Combiner c1;

    Combiner c2(l1,a1,l2);

    c2.init();

    CHECK(c2.check_init());
    CHECK(c2.hasindex(l1));
    CHECK(c2.hasindex(a1));
    CHECK(c2.hasindex(l2));
    CHECK_EQUAL(c2.right().m(),l1.m()*a1.m()*l2.m());

}

BOOST_AUTO_TEST_CASE(addLeft)
{
    Combiner c1;

    c1.addleft(l2);
    c1.addleft(a2);
    c1.addleft(l3);
    c1.addleft(l4);

    CHECK(c1.hasindex(l2));
    CHECK(c1.hasindex(l3));
    CHECK(c1.hasindex(a2));
    CHECK(c1.hasindex(l4));

    c1.init("cname");

    CHECK_EQUAL(c1.right().name(),"cname");
    CHECK_EQUAL(c1.right().m(),l2.m()*l3.m()*l4.m()*a2.m());

}

BOOST_AUTO_TEST_CASE(Product)
{

    ITensor A(a1,b3,l2,a4,l3);
    A.Randomize();

    Combiner c; 
    c.addleft(b3);
    c.addleft(l3);
    c.addleft(a1);

    c.init();
    Index r = c.right();

    ITensor cA = c * A;

    CHECK(cA.hasindex(l2));
    CHECK(cA.hasindex(a4));
    CHECK(cA.hasindex(r));
    CHECK(!cA.hasindex(b3));
    CHECK(!cA.hasindex(l3));
    CHECK(!cA.hasindex(a1));

    CHECK_CLOSE(A(l2(1),l3(1),b3(1)),cA(l2(1),r(1)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(1),b3(2)),cA(l2(1),r(2)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(1),b3(3)),cA(l2(1),r(3)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(2),b3(1)),cA(l2(1),r(4)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(2),b3(2)),cA(l2(1),r(5)),1E-10);
    CHECK_CLOSE(A(l2(1),l3(2),b3(3)),cA(l2(1),r(6)),1E-10);

    ITensor ucA = c * cA;

    CHECK(ucA.hasindex(a1));
    CHECK(ucA.hasindex(b3));
    CHECK(ucA.hasindex(l2));
    CHECK(ucA.hasindex(a4));
    CHECK(ucA.hasindex(l3));
    CHECK(!ucA.hasindex(r));

    for(int k2 = 1; k2 <= 2; ++k2)
    for(int k3 = 1; k3 <= 2; ++k3)
    for(int j3 = 1; j3 <= 2; ++j3)
    {
        CHECK_CLOSE(ucA(l2(k2),l3(k3),b3(j3)),A(l2(k2),l3(k3),b3(j3)),1E-10);
    }

    ITensor B(a1,b3,a3,l2,a4,l3);
    B.Randomize();

    Combiner c2; 
    c2.addleft(a3);
    c2.addleft(a1);

    ITensor cB = c2 * B;

    Index r2 = c2.right();

    CHECK(cB.hasindex(b3));
    CHECK(cB.hasindex(l2));
    CHECK(cB.hasindex(a4));
    CHECK(cB.hasindex(l3));
    CHECK(!cB.hasindex(a3));
    CHECK(!cB.hasindex(a1));

    for(int k2 = 1; k2 <= 2; ++k2)
    for(int k3 = 1; k3 <= 2; ++k3)
    for(int j3 = 1; j3 <= 2; ++j3)
    {
        CHECK_CLOSE(cB(b3(j3),l2(k2),l3(k3)),B(b3(j3),l2(k2),l3(k3)),1E-10);
    }

}

BOOST_AUTO_TEST_SUITE_END()
