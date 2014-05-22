#include "test.h"
#include "combiner.h"

TEST_CASE("CombinerTest")
{
    Index s1("s1",2,Site);
    Index s2("s2",2,Site);
    Index s3("s3",2,Site);
    Index s4("s4",2,Site);
    Index s1P(primed(s1));
    Index s2P(primed(s2));
    Index s3P(primed(s3));
    Index s4P(primed(s4));
    Index l1("l1",2);
    Index l2("l2",2);
    Index l3("l3",2);
    Index l4("l4",2);
    Index l5("l5",2);
    Index l6("l6",2);
    Index l7("l7",2);
    Index l8("l8",2);
    Index a1("a1");
    Index a2("a2");
    Index a3("a3");
    Index a4("a4");
    Index b2("b2",2);
    Index b3("b3",3);
    Index b4("b4",4);
    Index b5("b5",5);

SECTION("Constructors")
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

SECTION("addLeft")
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

SECTION("Product")
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

}
