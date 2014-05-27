#include "test.h"
#include "iqtensor.h"

using namespace itensor;
typedef IndexSet<IQIndex>
IQIndexSet;

TEST_CASE("IndexSetTest")
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
Index l30("Link3 Z0",1,Link);

IQIndex S1,S2,L1,L2,L3;

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

SECTION("Constructors")
    {
    shared_ptr<IQIndexSet> p1(new IQIndexSet(S1));
    CHECK_EQUAL(p1->index(1),S1);
    CHECK_CLOSE(p1->uniqueReal(),S1.uniqueReal(),1E-10);

    shared_ptr<IQIndexSet> p2(new IQIndexSet(S1,L1));
    CHECK_EQUAL(p2->index(1),S1);
    CHECK_EQUAL(p2->index(2),L1);
    const Real ur2 = S1.uniqueReal()
                   + L1.uniqueReal();
    CHECK_CLOSE(p2->uniqueReal(),ur2,1E-10);

    shared_ptr<IQIndexSet> p3(new IQIndexSet(S1,L1,S2));
    CHECK_EQUAL(p3->index(1),S1);
    CHECK_EQUAL(p3->index(2),L1);
    CHECK_EQUAL(p3->index(3),S2);
    const Real ur3 = S1.uniqueReal()
                   + L1.uniqueReal()
                   + S2.uniqueReal();
    CHECK_CLOSE(p3->uniqueReal(),ur3,1E-10);

    shared_ptr<IQIndexSet> p4(new IQIndexSet(S1,L1,S2,L2));
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

    shared_ptr<IQIndexSet> p5(new IQIndexSet(S1,L3,S2,L2));
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

SECTION("PrimeLevelMethods")
    {
    //
    // prime a specific IQIndex
    //
    shared_ptr<IQIndexSet> P = make_shared<IQIndexSet>(S1,prime(S1),S2,L2);

    P->prime(prime(S1),2);
    CHECK(P->index(1) == S1);
    CHECK(P->index(2) == prime(S1,3));
    CHECK(P->index(3) == S2);
    CHECK(P->index(4) == L2);
    }

SECTION("PrimeIndex")
    {
    shared_ptr<IQIndexSet> P = make_shared<IQIndexSet>(S1,prime(S2));

    P->prime(conj(S1));

    CHECK_EQUAL(P->index(1),prime(S1));
    //Even though the IQIndex passed to noprime had a different direction,
    //it still compares equal and the unprime IQIndex's arrow should be 
    //unchanged
    CHECK_EQUAL(P->index(1).dir(),S1.dir());
    }

SECTION("NoPrimeIndex")
    {
    shared_ptr<IQIndexSet> P = make_shared<IQIndexSet>(S1,prime(S2));

    P->noprime(conj(prime(S2)));

    CHECK_EQUAL(P->index(2),S2);
    //Even though the IQIndex passed to noprime had a different direction,
    //it still compares equal and the unprime IQIndex's arrow should be 
    //unchanged
    CHECK_EQUAL(P->index(2).dir(),S2.dir());
    }

SECTION("NoPrimeType")
    {
    IQIndexSet I1(S1,prime(S2),prime(L1),prime(L2)),
               I2(S1,prime(S2),L1,prime(L1));

    I1.noprime(Link);

    CHECK_THROWS_AS(I2.noprime(Link),ITError);
    }

SECTION("AddIndex")
    {
    shared_ptr<IQIndexSet> P = make_shared<IQIndexSet>();

    P->addindex(S1);
    P->addindex(prime(S1));
    P->addindex(L2);

    CHECK_EQUAL(P->index(1),S1);
    CHECK_EQUAL(P->index(2),prime(S1));
    CHECK_EQUAL(P->index(3),L2);
    CHECK_EQUAL(P->r(),3);
    }

SECTION("Contraction")
    {
    IQIndexSet is1(S1,L1,L3),
               is2(L1,L2,L3,S2);

    IQIndexSet res = is1 * is2;

    CHECK(hasindex(res,S1));
    CHECK(!hasindex(res,L1));
    CHECK(!hasindex(res,L3));
    CHECK(hasindex(res,L2));
    CHECK(hasindex(res,S2));

    IndexSet<Index> is3(s1u,l1u,l2uu,s1d),
                    is4(s1d,s2u,l2dd,l1u);

    IndexSet<Index> res2 = is3 * is4;

    CHECK(hasindex(res2,s1u));
    CHECK(!hasindex(res2,l1u));
    CHECK(hasindex(res2,l2uu));
    CHECK(!hasindex(res2,s1d));
    CHECK(hasindex(res2,s2u));
    CHECK(hasindex(res2,l2dd));
    }

}
