#include "test.h"
#include "iqcombiner.h"

using namespace itensor;
using namespace std;

TEST_CASE("IQCombinerTest")
{

Index s1u("Site1 Up",1,Site);
Index s1d("Site1 Dn",1,Site);
Index s2u("Site2 Up",1,Site);
Index s2d("Site2 Dn",1,Site);
Index l1u("Link1 Up",2,Link);
Index l10("Link1 Z0",2,Link);
Index l1d("Link1 Dn",2,Link);
Index l2uu("Link2 UU",2,Link);
Index l20("Link2 Z0",2,Link);
Index l2dd("Link2 DD",2,Link);

IQIndex S1,S2,L1,L2;

IQTensor phi;

S1 = IQIndex("S1",
             s1u,QN(+1),
             s1d,QN(-1),Out);
S2 = IQIndex("S2",
             s2u,QN(+1),
             s2d,QN(-1),Out);
L1 = IQIndex("L1",
             l1u,QN(+1),
             l10,QN( 0),
             l1d,QN(-1),
             Out);
L2 = IQIndex("L2",
             l2uu,QN(+2),
             l20,QN( 0),
             l2dd,QN(-2),
             Out);

phi = IQTensor(S1,S2,L2);
    {
    ITensor uu(s1u,s2u,l2dd);
    uu.randomize();
    phi += uu;

    ITensor ud(s1u,s2d,l20);
    ud.randomize();
    phi += ud;

    ITensor du(s1d,s2u,l20);
    du.randomize();
    phi += du;
    }

SECTION("Constructors")
    {
    IQCombiner c1;

    IQCombiner c2(L1,S1,L2);

    c2.init();

    CHECK(c2.isInit());
    CHECK(hasindex(c2,L1));
    CHECK(hasindex(c2,S1));
    CHECK(hasindex(c2,L2));
    CHECK_EQUAL(c2.right().m(),L1.m()*S1.m()*L2.m());
    }

SECTION("addLeft")
    {
    IQCombiner c1,
               c1s;

    c1.addleft(L2);
    c1.addleft(S2);
    c1.addleft(L1);
    c1.addleft(S1);
    c1.init("cname");

    c1s.addleft(L2);
    c1s.addleft(S2);
    c1s.addleft(L1);
    c1s.addleft(S1);
    c1s.init("cname",Site);

    CHECK(hasindex(c1,L2));
    CHECK(hasindex(c1,S2));
    CHECK(hasindex(c1,L1));
    CHECK(hasindex(c1,S1));

    //CHECK_EQUAL(c1.right().name(),"cname");
    //CHECK_EQUAL(c1s.right().name(),"cname");
    //CHECK_EQUAL(c1.right().type(),Link);
    //CHECK_EQUAL(c1s.right().type(),Site);
    CHECK_EQUAL(c1.right().m(),L2.m()*S2.m()*L1.m()*S1.m());
    CHECK_EQUAL(c1s.right().m(),L2.m()*S2.m()*L1.m()*S1.m());

    }

SECTION("Product")
    {
    CHECK_EQUAL(L2.dir(),Out);

    IQCombiner c;
    c.addleft(S2);
    c.addleft(L2);
    c.init();

    IQTensor cphi = c * phi;

    CHECK(hasindex(cphi,S1));
    CHECK(hasindex(cphi,c.right()));

    IQIndex r = c.right();

    //Now uncombine - to do so
    //must use the conjugate IQCombiner

    IQCombiner cc(c);
    cc.conj();
    CHECK(c.right().dir() != cc.right().dir());

    IQTensor ucphi = cc * cphi;

    CHECK(hasindex(ucphi,S1));
    CHECK(hasindex(ucphi,S2));
    CHECK(hasindex(ucphi,L2));

    IQTensor diff = phi - ucphi;
    CHECK(diff.norm() < 1E-12);
    }

SECTION("Primes")
    {
    IQCombiner c;
    c.addleft(S2);
    c.addleft(L2);
    c.init();

    IQTensor pphi = prime(phi);

    IQTensor cpphi = prime(c) * pphi;

    CHECK(hasindex(cpphi,prime(S1)));
    CHECK(hasindex(cpphi,prime(c.right())));

    //Check that using a prime combiner gives
    //same result as regular combiner, then
    //priming
    IQTensor cphi = c * phi;
    IQTensor diff = prime(cphi) - cpphi;
    CHECK(diff.norm() < 1E-10);

    }

SECTION("CondenseProduct")
    {
    CHECK_EQUAL(L2.dir(),Out);

    IQTensor phi(S1,S2,L2);

    ITensor ud(s1u,s2d,l20);
    ud.randomize();
    phi += ud;

    ITensor du(s1d,s2u,l20);
    du.randomize();
    phi += du;

    const QN Zero;
    CHECK(div(phi) == Zero);

    IQCombiner c;
    c.doCondense(true);
    CHECK(c.doCondense());

    c.addleft(S2);
    c.addleft(L2);
    c.init();

    IQTensor cphi = c * phi;

    CHECK(hasindex(cphi,S1));
    CHECK(hasindex(cphi,c.right()));

    IQCombiner cc(c);
    cc.conj();
    CHECK(c.right().dir() != cc.right().dir());

    IQTensor ucphi = conj(c) * cphi;

    CHECK(hasindex(ucphi,S1));
    CHECK(hasindex(ucphi,S2));
    CHECK(hasindex(ucphi,L2));

    IQTensor diff = phi - ucphi;
    CHECK(diff.norm() < 1E-12);
    }
}

