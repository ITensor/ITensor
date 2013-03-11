#include "test.h"
#include "iqcombiner.h"
#include <boost/test/unit_test.hpp>

using namespace std;

struct IQCombinerDefaults
    {
    const Index
    s1u,s1d,s2u,s2d,
    l1u,l10,l1d,
    l2uu,l2u,l20,l2d,l2dd;

    IQIndex S1,S2,L1,L2;

    IQTensor phi;

    IQCombinerDefaults() :
        s1u(Index("Site1 Up",1,Site)),
        s1d(Index("Site1 Dn",1,Site)),
        s2u(Index("Site2 Up",1,Site)),
        s2d(Index("Site2 Dn",1,Site)),
        l1u(Index("Link1 Up",2,Link)),
        l10(Index("Link1 Z0",2,Link)),
        l1d(Index("Link1 Dn",2,Link)),
        l2uu(Index("Link2 UU",2,Link)),
        l2u(Index("Link2 Up",2,Link)),
        l20(Index("Link2 Z0",2,Link)),
        l2d(Index("Link2 Dn",2,Link)),
        l2dd(Index("Link2 DD",2,Link))
        {
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
                     l2u,QN(+1),
                     l20,QN( 0),
                     l2d,QN(-1),
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
        }

    ~IQCombinerDefaults() { }

    };

BOOST_FIXTURE_TEST_SUITE(IQCombinerTest,IQCombinerDefaults)

TEST(Constructors)
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

TEST(addLeft)
    {
    IQCombiner c1;

    c1.addleft(L2);
    c1.addleft(S2);
    c1.addleft(L1);
    c1.addleft(S1);

    CHECK(hasindex(c1,L2));
    CHECK(hasindex(c1,S2));
    CHECK(hasindex(c1,L1));
    CHECK(hasindex(c1,S1));

    IQCombiner c1s(c1);

    c1.init("cname");
    c1s.init("cname",Site);

    CHECK_EQUAL(c1.right().name(),"cname");
    CHECK_EQUAL(c1s.right().name(),"cname");
    CHECK_EQUAL(c1.right().type(),Link);
    CHECK_EQUAL(c1s.right().type(),Site);
    CHECK_EQUAL(c1.right().m(),L2.m()*S2.m()*L1.m()*S1.m());
    CHECK_EQUAL(c1s.right().m(),L2.m()*S2.m()*L1.m()*S1.m());

    }

TEST(Product)
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
    CHECK_CLOSE(diff.norm(),0,1E-10);
    }

TEST(Primes)
    {
    IQCombiner c;
    c.addleft(S2);
    c.addleft(L2);
    c.init();

    IQTensor pphi = primed(phi);

    IQTensor cpphi = primed(c) * pphi;

    CHECK(hasindex(cpphi,primed(S1)));
    CHECK(hasindex(cpphi,primed(c.right())));

    //Check that using a primed combiner gives
    //same result as regular combiner, then
    //priming
    IQTensor cphi = c * phi;
    IQTensor diff = primed(cphi) - cpphi;
    CHECK(diff.norm() < 1E-10);

    }

TEST(CondenseProduct)
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

BOOST_AUTO_TEST_SUITE_END()
