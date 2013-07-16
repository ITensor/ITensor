#include "test.h"
#include "iqtsparse.h"
#include <boost/test/unit_test.hpp>

struct IQTSparseDefaults
    {
    const Index
    s1u,s1d,s2u,s2d,
    l1u,l10,l1d,
    l2uu,l20,l2dd;

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
    l20(Index("Link2 Z0",2,Link)),
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
                     l20,QN( 0),
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

        A = IQTensor(L1(1),S1(1),L2(1),S2(1));
        A.randomize();

        B = IQTensor(L1(l1u.m()+1),L2(l2uu.m()+1));
        B.randomize();

        C = IQTensor(conj(L1)(l1u.m()+1),primed(L1)(l1u.m()+1));
        C.randomize();

        D = IQTensor(conj(L1)(l1u.m()-1),S1(2),primed(L1)(l1u.m()+1),primed(L1,2)(l1u.m()+1));
        D.randomize();
        }

    void
    test1()
        {
        IQTSparse B(L1,L2);

        Print(hasindex(B,L1));
        Print(hasindex(B,L2));
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

    DD += ITSparse(l2uu,l1d,2);
    }

TEST(ContractingProduct)
    {
    }


BOOST_AUTO_TEST_SUITE_END()
