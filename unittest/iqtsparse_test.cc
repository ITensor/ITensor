#include "test.h"
#include "iqtsparse.h"
#include <boost/test/unit_test.hpp>

struct IQTSparseDefaults
    {
    const Index
    s1u,s1d,s2u,s2d,
    l1u,l10,l1d,
    l2uu,l2u,l20,l2d,l2dd;

    IQIndex S1,S2,L1,L2;

    IQTensor phi,A,B,C,D;

    IQTSparseDefaults() :
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

BOOST_FIXTURE_TEST_SUITE(IQTSparseTest,IQTSparseDefaults)

TEST(Constructors)
    {
    IQTSparse B(L1,L2);

    CHECK(hasindex(B,L1));
    CHECK(hasindex(B,L2));
    }

TEST(Addition)
    {
    IQTSparse DD(L1,L2);

    DD += ITSparse(l2u,l1d,2);
    }

TEST(ContractingProduct)
    {
    }


BOOST_AUTO_TEST_SUITE_END()
