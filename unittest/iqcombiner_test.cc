#include "test.h"
#include "iqcombiner.h"
#include <boost/test/unit_test.hpp>

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
            uu.Randomize();
            phi += uu;

            ITensor ud(s1u,s2d,l20);
            ud.Randomize();
            phi += ud;

            ITensor du(s1d,s2u,l20);
            du.Randomize();
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
    CHECK(c2.hasindex(L1));
    CHECK(c2.hasindex(S1));
    CHECK(c2.hasindex(L2));
    CHECK_EQUAL(c2.right().m(),L1.m()*S1.m()*L2.m());
    }

TEST(addLeft)
    {
    IQCombiner c1;

    c1.addleft(L2);
    c1.addleft(S2);
    c1.addleft(L1);
    c1.addleft(S1);

    CHECK(c1.hasindex(L2));
    CHECK(c1.hasindex(S2));
    CHECK(c1.hasindex(L1));
    CHECK(c1.hasindex(S1));

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

    CHECK(cphi.hasindex(S1));
    CHECK(cphi.hasindex(c.right()));

    IQIndex r = c.right();

    //Now uncombine - to do so
    //must use the conjugate IQCombiner

    IQCombiner cc(c);
    cc.conj();
    CHECK(c.right().dir() != cc.right().dir());

    IQTensor ucphi = cc * cphi;

    CHECK(ucphi.hasindex(S1));
    CHECK(ucphi.hasindex(S2));
    CHECK(ucphi.hasindex(L2));

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

    CHECK(cpphi.hasindex(primed(S1)));
    CHECK(cpphi.hasindex(primed(c.right())));

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

    ITensor uu(s1u,s2u,l2dd);
    uu.Randomize();
    phi += uu;

    ITensor ud(s1u,s2d,l20);
    ud.Randomize();
    phi += ud;

    ITensor du(s1d,s2u,l20);
    du.Randomize();
    phi += du;

    checkDiv(phi);

    IQCombiner c;
    c.doCondense(true);
    CHECK(c.doCondense());

    c.addleft(S2);
    c.addleft(L2);
    c.init();

    IQTensor cphi = c * phi;

    CHECK(cphi.hasindex(S1));
    CHECK(cphi.hasindex(c.right()));

    IQCombiner cc(c);
    cc.conj();
    CHECK(c.right().dir() != cc.right().dir());

    IQTensor ucphi = conj(c) * cphi;

    CHECK(ucphi.hasindex(S1));
    CHECK(ucphi.hasindex(S2));
    CHECK(ucphi.hasindex(L2));

    IQTensor diff = phi - ucphi;
    CHECK(diff.norm() < 1E-12);
    }

BOOST_AUTO_TEST_SUITE_END()
