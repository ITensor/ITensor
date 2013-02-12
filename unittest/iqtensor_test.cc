#include "test.h"
#include "iqtensor.h"
#include <boost/test/unit_test.hpp>

using namespace std;

struct IQTensorDefaults
    {
    const Index
    s1u,s1d,s2u,s2d,
    l1u,l10,l1d,
    l2uu,l2u,l20,l2d,l2dd;

    IQIndex S1,S2,L1,L2;

    IQTensor phi,A,B,C,D;

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
        uu.randomize();
        phi += uu;

        ITensor ud(s1u,s2d,l20);
        ud.randomize();
        phi += ud;

        ITensor du(s1d,s2u,l20);
        du.randomize();
        phi += du;
        }

        A = IQTensor(L1,S1,L2,S2);
        for(int n1 = 1; n1 <= L1.nindex(); ++n1)
        for(int n2 = 1; n2 <= L2.nindex(); ++n2)
        for(int p1 = 1; p1 <= S1.nindex(); ++p1)
        for(int p2 = 1; p2 <= S2.nindex(); ++p2)
            {
            ITensor T(L1.index(n1),L2.index(n2),S1.index(p1),S2.index(p2));
            T.randomize();
            A += T;
            }

        B = IQTensor(L1,L2);
        for(int n1 = 1; n1 <= L1.nindex(); ++n1)
        for(int n2 = 1; n2 <= L2.nindex(); ++n2)
            {
            ITensor T(L1.index(n1),L2.index(n2));
            T.randomize();
            B += T;
            }

        C = IQTensor(conj(L1),primed(L1));
        for(int n1 = 1; n1 <= L1.nindex(); ++n1)
            {
            Matrix U(L1.index(n1).m(),L1.index(n1).m());
            U.Randomize();
            U += U.t();
            ITensor T(L1.index(n1),primed(L1).index(n1),U);
            C += T;
            }

        D = IQTensor(conj(L1),S1,primed(L1),primed(L1,2));
        for(int n1 = 1; n1 <= L1.nindex(); ++n1)
        for(int n2 = 1; n2 <= S1.nindex(); ++n2)
        for(int n3 = 1; n3 <= L1.nindex(); ++n3)
        for(int n4 = 1; n4 <= S1.nindex(); ++n4)
            {
            ITensor T(L1.index(n1),S1.index(n2),primed(L1).index(n3),primed(L1,2).index(n4));
            T.randomize();
            D += T;
            }
        }

    };

BOOST_FIXTURE_TEST_SUITE(IQTensorTest,IQTensorDefaults)

TEST(Null)
    {
    IQTensor t1;

    CHECK(t1.isNull());
    }

TEST(Constructors)
    {
    Real f = ran1();
    IQTensor rZ(f);

    CHECK_EQUAL(rZ.r(),0);
    CHECK_CLOSE(rZ.norm(),f,1E-10);
    }

TEST(NonContractProd)
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

TEST(ITensorConversion)
    {

    ITensor itphi = phi.toITensor();

    for(int k1 = 1; k1 <= S1.m(); ++k1)
    for(int k2 = 1; k2 <= S2.m(); ++k2)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
        CHECK_CLOSE(phi(S1(k1),S2(k2),L2(j2)),itphi(Index(S1)(k1),Index(S2)(k2),Index(L2)(j2)),1E-5);

    ITensor itA = A.toITensor();

    for(int k1 = 1; k1 <= S1.m(); ++k1)
    for(int k2 = 1; k2 <= S2.m(); ++k2)
    for(int j1 = 1; j1 <= L1.m(); ++j1)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
        CHECK_CLOSE(A(S1(k1),S2(k2),L1(j1),L2(j2)),
                    itA(Index(S1)(k1),Index(S2)(k2),Index(L1)(j1),Index(L2)(j2)),
                    1E-5);


    }

TEST(SymmetricDiag11)
    {
    IQTensor D,U;
    IQIndex mid;
    int mink,maxk;
    C *= -2;
    C.symmetricDiag11(L1,D,U,mid,mink,maxk);

    IQTensor UD(U);
    UD.prime(L1);
    UD /= D;
    ITensor diff = (conj(UD)*U - C).toITensor();
    CHECK(diff.norm() < 1E-10);

    IQTensor set1(conj(mid));
    set1(mid(mink)) = 1;
    U *= set1;
    CHECK_CLOSE(D(mid(mink)),Dot(primed(U),C*U),1E-10);

    }

TEST(TieIndices)
    {
    IQTensor D1 = tieIndices(D,L1,L1.primed(),L1);

    for(int k1 = 1; k1 <= L1.m(); ++k1)
    for(int k2 = 1; k2 <= L1.m(); ++k2)
    for(int k3 = 1; k3 <= S1.m(); ++k3)
        {
        CHECK_CLOSE(D1(L1(k1),primed(L1,2)(k2),S1(k3)),D(L1(k1),primed(L1)(k1),primed(L1,2)(k2),S1(k3)),1E-10);
        }
    }

TEST(ToReal)
    {
    Real f = ran1();
    IQTensor T(f);
    CHECK_CLOSE(T.toReal(),f,1E-5);

    //Default constructed IQTensor should
    //have a toReal value of zero
    IQTensor Z;
    CHECK_CLOSE(Z.toReal(),0,1E-5);
    }

TEST(Trace)
    {

    Real f = -ran1();
    D *= f;

    IQTensor Dt = trace(conj(L1),primed(L1,2),D);

    for(int j2 = 1; j2 <= S1.m(); ++j2)
    for(int j1 = 1; j1 <= L1.m(); ++j1)
        {
        Real val = 0;
        for(int k1 = 1; k1 <= L1.m(); ++k1)
            {
            val += D(conj(L1)(k1),S1(j2),primed(L1)(j1),primed(L1,2)(k1));
            }
        CHECK_CLOSE(val,Dt(S1(j2),primed(L1)(j1)),1E-10);
        }
    }

TEST(MapElems)
    {
    IQTensor B1(B);

    Functor f;
    B1.mapElems(f);

    for(int j1 = 1; j1 <= L1.m(); ++j1)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
        {
        CHECK_CLOSE( f( B(L1(j1),L2(j2)) ), 
                    B1(L1(j1),L2(j2)),1E-3);
        }
    }

TEST(Randomize)
    {
    IQTensor T(S1,S2,L2);
    T += ITensor(s1u,s2u,l2dd);
    T.randomize();
    //PrintDat(T);
    }

BOOST_AUTO_TEST_SUITE_END()
