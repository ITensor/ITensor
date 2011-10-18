#include "test.h"
#include "iqtensor.h"
#include <boost/test/unit_test.hpp>

struct IQTensorDefaults
{
    const Index
    s1u,s1d,s2u,s2d,
    l1u,l10,l1d,
    l2uu,l2u,l20,l2d,l2dd;

    IQIndex S1,S2,L1,L2;

    IQTensor phi,A,B;

    IQTensorDefaults() :
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
        {
        phi = IQTensor(S1,S2,L2);

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

        A = IQTensor(L1,S1,L2,S2);
        for(int n1 = 1; n1 <= L1.nindex(); ++n1)
        for(int n2 = 1; n2 <= L2.nindex(); ++n2)
        for(int p1 = 1; p1 <= S1.nindex(); ++p1)
        for(int p2 = 1; p2 <= S2.nindex(); ++p2)
            {
            ITensor T(L1.index(n1),L2.index(n2),S1.index(p1),S2.index(p2));
            T.Randomize();
            A += T;
            }

        B = IQTensor(L1,L2);
        for(int n1 = 1; n1 <= L1.nindex(); ++n1)
        for(int n2 = 1; n2 <= L2.nindex(); ++n2)
            {
            ITensor T(L1.index(n1),L2.index(n2));
            T.Randomize();
            B += T;
            }
    }

};

BOOST_FIXTURE_TEST_SUITE(IQTensorTest,IQTensorDefaults)

BOOST_AUTO_TEST_CASE(Null)
    {
    IQTensor t1;

    CHECK(t1.is_null());
    }

BOOST_AUTO_TEST_CASE(Constructors)
    {
    Real f = ran1();
    IQTensor rZ(f);

    CHECK_EQUAL(rZ.r(),0);
    CHECK_CLOSE(rZ.norm(),f,1E-10);
    }

BOOST_AUTO_TEST_CASE(NonContractProd)
    {

    IQTensor res = A / B;

    for(int j1 = 1; j1 <= L1.m(); ++j1)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
    for(int k1 = 1; k1 <= S1.m(); ++k1)
    for(int k2 = 1; k2 <= S2.m(); ++k2)
        {
        CHECK_CLOSE(res(L1(j1),L2(j2),S1(k1),S2(k2)),
                    A(L1(j1),S1(k1),L2(j2),S2(k2))*B(L1(j1),L2(j2)),1E-5);
        }

    }

BOOST_AUTO_TEST_CASE(ITensorConversion)
    {

    ITensor itphi = phi;

    for(int k1 = 1; k1 <= S1.m(); ++k1)
    for(int k2 = 1; k2 <= S2.m(); ++k2)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
        CHECK_CLOSE(phi(S1(k1),S2(k2),L2(j2)),itphi(Index(S1)(k1),Index(S2)(k2),Index(L2)(j2)),1E-5);

    ITensor itA = A;

    for(int k1 = 1; k1 <= S1.m(); ++k1)
    for(int k2 = 1; k2 <= S2.m(); ++k2)
    for(int j1 = 1; j1 <= L1.m(); ++j1)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
        CHECK_CLOSE(A(S1(k1),S2(k2),L1(j1),L2(j2)),
                    itA(Index(S1)(k1),Index(S2)(k2),Index(L1)(j1),Index(L2)(j2)),
                    1E-5);


    }

BOOST_AUTO_TEST_SUITE_END()
