#include "test.h"
#include "iqtensor.h"
#include <boost/test/unit_test.hpp>
#include <boost/intrusive_ptr.hpp>

using namespace std;
using namespace boost;

typedef IndexSet<IQIndex>
IQIndexSet;

struct IndexSetDefaults
    {
    const Index
    s1u,s1d,s2u,s2d,
    l1u,l10,l1d,
    l2uu,l20,l2dd,
    l30;

    IQIndex S1,S2,L1,L2,L3;

    IndexSetDefaults() :
    s1u(Index("Site1 Up",1,Site)),
    s1d(Index("Site1 Dn",1,Site)),
    s2u(Index("Site2 Up",1,Site)),
    s2d(Index("Site2 Dn",1,Site)),
    l1u(Index("Link1 Up",2,Link)),
    l10(Index("Link1 Z0",2,Link)),
    l1d(Index("Link1 Dn",2,Link)),
    l2uu(Index("Link2 UU",2,Link)),
    l20(Index("Link2 Z0",2,Link)),
    l2dd(Index("Link2 DD",2,Link)),
    l30(Index("Link3 Z0",1,Link))
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
        L3 = IQIndex("L3",
                     l30,QN(0),
                     In);
        }
    };

BOOST_FIXTURE_TEST_SUITE(IndexSetTest,IndexSetDefaults)

TEST(Constructors)
    {
    boost::shared_ptr<IQIndexSet> p1(new IQIndexSet(S1));
    CHECK_EQUAL(p1->index(1),S1);
    CHECK_CLOSE(p1->uniqueReal(),S1.uniqueReal(),1E-10);

    boost::shared_ptr<IQIndexSet> p2(new IQIndexSet(S1,L1));
    CHECK_EQUAL(p2->index(1),S1);
    CHECK_EQUAL(p2->index(2),L1);
    const Real ur2 = S1.uniqueReal()
                   + L1.uniqueReal();
    CHECK_CLOSE(p2->uniqueReal(),ur2,1E-10);

    boost::shared_ptr<IQIndexSet> p3(new IQIndexSet(S1,L1,S2));
    CHECK_EQUAL(p3->index(1),S1);
    CHECK_EQUAL(p3->index(2),L1);
    CHECK_EQUAL(p3->index(3),S2);
    const Real ur3 = S1.uniqueReal()
                   + L1.uniqueReal()
                   + S2.uniqueReal();
    CHECK_CLOSE(p3->uniqueReal(),ur3,1E-10);

    boost::shared_ptr<IQIndexSet> p4(new IQIndexSet(S1,L1,S2,L2));
    CHECK_EQUAL(p4->index(1), S1);
    CHECK_EQUAL(p4->index(2), L1);
    CHECK_EQUAL(p4->index(3), S2);
    CHECK_EQUAL(p4->index(4), L2);
    const Real ur4 = S1.uniqueReal()
                   + L1.uniqueReal()
                   + S2.uniqueReal()
                   + L2.uniqueReal();
    CHECK_CLOSE(p4->uniqueReal(),ur4,1E-10);

    //Check that m==1 indices get sorted to the back

    CHECK_EQUAL(L3.m(),1);

    boost::shared_ptr<IQIndexSet> p5(new IQIndexSet(S1,L3,S2,L2));
    CHECK_EQUAL(p5->index(1),S1);
    CHECK_EQUAL(p5->index(2),S2);
    CHECK_EQUAL(p5->index(3),L2);
    CHECK_EQUAL(p5->index(4),L3);
    const Real ur5 = S1.uniqueReal()
                   + S2.uniqueReal()
                   + L2.uniqueReal()
                   + L3.uniqueReal();
    CHECK_CLOSE(p5->uniqueReal(),ur5,1E-10);
    }

TEST(PrimeLevelMethods)
    {
    //
    // prime a specific IQIndex
    //
    boost::shared_ptr<IQIndexSet> P = boost::make_shared<IQIndexSet>(S1,primed(S1),S2,L2);

    P->prime(primed(S1),2);
    CHECK(P->index(1) == S1);
    CHECK(P->index(2) == primed(S1,3));
    CHECK(P->index(3) == S2);
    CHECK(P->index(4) == L2);
    }

TEST(PrimeIndex)
    {
    boost::shared_ptr<IQIndexSet> P = boost::make_shared<IQIndexSet>(S1,primed(S2));

    P->prime(conj(S1));

    CHECK_EQUAL(P->index(1),primed(S1));
    //Even though the IQIndex passed to noprime had a different direction,
    //it still compares equal and the unprimed IQIndex's arrow should be 
    //unchanged
    CHECK_EQUAL(P->index(1).dir(),S1.dir());
    }

TEST(NoPrimeIndex)
    {
    boost::shared_ptr<IQIndexSet> P = boost::make_shared<IQIndexSet>(S1,primed(S2));

    P->noprime(conj(primed(S2)));

    CHECK_EQUAL(P->index(2),S2);
    //Even though the IQIndex passed to noprime had a different direction,
    //it still compares equal and the unprimed IQIndex's arrow should be 
    //unchanged
    CHECK_EQUAL(P->index(2).dir(),S2.dir());
    }

TEST(NoPrimeType)
    {
    IQIndexSet I1(S1,primed(S2),primed(L1),primed(L2)),
               I2(S1,primed(S2),L1,primed(L1));

    I1.noprime(Link);

    CHECK_THROW(I2.noprime(Link),ITError);
    }

TEST(AddIndex)
    {
    boost::shared_ptr<IQIndexSet> P = boost::make_shared<IQIndexSet>();

    P->addindex(S1);
    P->addindex(primed(S1));
    P->addindex(L2);

    CHECK_EQUAL(P->index(1),S1);
    CHECK_EQUAL(P->index(2),primed(S1));
    CHECK_EQUAL(P->index(3),L2);
    CHECK_EQUAL(P->r(),3);
    }

BOOST_AUTO_TEST_SUITE_END()
