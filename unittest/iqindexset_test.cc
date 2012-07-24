#include "test.h"
#include "iqtensor.h"
#include <boost/test/unit_test.hpp>
#include <boost/intrusive_ptr.hpp>

using namespace std;
using namespace boost;

struct IQIndexSetDefaults
    {
    const Index
    s1u,s1d,s2u,s2d,
    l1u,l10,l1d,
    l2uu,l2u,l20,l2d,l2dd;

    IQIndex S1,S2,L1,L2;


    IQIndexSetDefaults() :
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
        }

    };

BOOST_FIXTURE_TEST_SUITE(IQIndexSetTest,IQIndexSetDefaults)

TEST(Constructors)
    {
    intrusive_ptr<IQIndexSet> p1 = new IQIndexSet(S1);
    CHECK(p1->index(1) == S1);

    intrusive_ptr<IQIndexSet> p2 = new IQIndexSet(S1,L1);
    CHECK(p2->index(1) == S1);
    CHECK(p2->index(2) == L1);

    intrusive_ptr<IQIndexSet> p3 = new IQIndexSet(S1,L1,S2);
    CHECK(p3->index(1) == S1);
    CHECK(p3->index(2) == L1);
    CHECK(p3->index(3) == S2);

    intrusive_ptr<IQIndexSet> p4 = new IQIndexSet(S1,L1,S2,L2);
    CHECK(p4->index(1) == S1);
    CHECK(p4->index(2) == L1);
    CHECK(p4->index(3) == S2);
    CHECK(p4->index(4) == L2);
    }

TEST(PrimeLevelMethods)
    {
    //
    // indIncPrime
    //
    intrusive_ptr<IQIndexSet> P = new IQIndexSet(S1,primed(S1),S2,L2);

    P->indIncPrime(primed(S1),2);
    CHECK(P->index(1) == S1);
    CHECK(P->index(2) == primed(S1,3));
    CHECK(P->index(3) == S2);
    CHECK(P->index(4) == L2);

    //
    // indIncAllPrime
    //
    P = new IQIndexSet(S1,primed(S1),S2,L2);

    P->indIncAllPrime(S1,2);
    CHECK(P->index(1) == primed(S1,2));
    CHECK(P->index(2) == primed(S1,3));
    CHECK(P->index(3) == S2);
    CHECK(P->index(4) == L2);
    }

BOOST_AUTO_TEST_SUITE_END()
