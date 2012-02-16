#include "test.h"
#include "itsparse.h"
#include <boost/test/unit_test.hpp>

struct ITSparseDefaults
    {
    Index s1,s2,s3,s4,
          s1P,s2P,s3P,s4P,
          l1,l2,l3,l4,l5,l6,l7,l8,
          a1,a2,a3,a4,
          b2,b3,b4,b5;

    ITensor A,B,X,Z;

    std::vector<Index> mixed_inds, reordered_mixed_inds;

    int mixed_inds_dim;

    ITSparseDefaults()
        :
        s1(Index("s1",2,Site)),
        s2(Index("s2",2,Site)),
        s3(Index("s3",2,Site)),
        s4(Index("s4",2,Site)),
        s1P(s1.primed()),
        s2P(s2.primed()),
        s3P(s3.primed()),
        s4P(s4.primed()),
        l1(Index("l1",2)),
        l2(Index("l2",2)),
        l3(Index("l3",2)),
        l4(Index("l4",2)),
        l5(Index("l5",2)),
        l6(Index("l6",2)),
        l7(Index("l7",2)),
        l8(Index("l8",2)),
        a1(Index("a1")),
        a2(Index("a2")),
        a3(Index("a3")),
        a4(Index("a4")),
        b2(Index("b2",2)),
        b3(Index("b3",3)),
        b4(Index("b4",4)),
        b5(Index("b5",5)),
        mixed_inds(6),
        reordered_mixed_inds(6),
        mixed_inds_dim(1)
    {
        {
        Matrix M(s1.m(),s2.m());
        M(1,1) = 11; M(1,2) = 12;
        M(2,1) = 21; M(2,2) = 22;
        A = ITensor(s1,s2,M);
        }

        {
        Matrix M(s1.m(),s2.m());
        M(1,1) = 110; M(1,2) = 120;
        M(2,1) = 210; M(2,2) = 220;
        B = ITensor(s1,s2,M);
        }

        {
        Matrix M(s1.m(),s2.m());
        M(1,1) = 0; M(1,2) = 1;
        M(2,1) = 1; M(2,2) = 0;
        X = ITensor(s1,s2,M);
        }
        
        {
        Matrix M(s1.m(),s2.m());
        M(1,1) = 1; M(1,2) =  0;
        M(2,1) = 0; M(2,2) = -1;
        Z = ITensor(s1,s2,M);
        }

        mixed_inds[0] = a2;
        mixed_inds[1] = b3;
        mixed_inds[2] = l1;
        mixed_inds[3] = l2;
        mixed_inds[4] = a4;
        mixed_inds[5] = l4;

        Foreach(const Index& I, mixed_inds)
        { mixed_inds_dim *= I.m(); }

        reordered_mixed_inds[0] = a2;
        reordered_mixed_inds[1] = l1;
        reordered_mixed_inds[2] = b3;
        reordered_mixed_inds[3] = a4;
        reordered_mixed_inds[4] = l4;
        reordered_mixed_inds[5] = l2;
    }

    ~ITSparseDefaults() { }

    }; //struct ITSparseDefaults

BOOST_FIXTURE_TEST_SUITE(ITSparseTest,ITSparseDefaults)

TEST(Constructors)
    {
    ITSparse B(b2,b3,2);

    CHECK_CLOSE(sqrt(2)*2,B.norm(),1E-5);
    CHECK(B.hasindex(b2));
    CHECK(B.hasindex(b3));
    CHECK_EQUAL(b2.m(),B.diagSize());

    Vector diag(b3.m());
    diag.Randomize();

    ITSparse D(b3,primed(b3),diag);
    D *= -1;
    CHECK(D.hasindex(b3));
    CHECK(D.hasindex(primed(b3)));
    CHECK_EQUAL(diag.Length(),D.diagSize());
    CHECK_CLOSE(Norm(diag),D.norm(),1E-5);
    }

TEST(ContractingProduct)
    {
    ITSparse D(b3,b4,1);

    ITensor T(a1,b2,b4,s1);
    T.Randomize();

    ITensor R = D * T;

    CHECK(!R.hasindex(b4));
    CHECK(R.hasindex(b3));
    CHECK(R.hasindex(b2));
    CHECK(R.hasindex(s1));
    CHECK(R.hasindex(a1));

    for(int j2 = 1; j2 <= b2.m(); ++j2)
    for(int j3 = 1; j3 <= b3.m(); ++j3)
    for(int k1 = 1; k1 <= s1.m(); ++k1)
        {
        CHECK_CLOSE(T(b2(j2),b4(j3),s1(k1)),R(b3(j3),b2(j2),s1(k1)),1E-5);
        }

    Vector diag(min(b3.m(),b4.m()));
    diag.Randomize();

    ITSparse D2(b3,b4,diag);

    ITensor T2(a1,b2,b4,s1);
    T2.Randomize();

    ITensor R2 = D2 * T2;

    CHECK(!R2.hasindex(b4));
    CHECK(R2.hasindex(b3));
    CHECK(R2.hasindex(b2));
    CHECK(R2.hasindex(s1));
    CHECK(R2.hasindex(a1));

    for(int j2 = 1; j2 <= b2.m(); ++j2)
    for(int j3 = 1; j3 <= b3.m(); ++j3)
    for(int k1 = 1; k1 <= s1.m(); ++k1)
        {
        CHECK_CLOSE(diag(j3) * T2(b2(j2),b4(j3),s1(k1)),R2(b3(j3),b2(j2),s1(k1)),1E-5);
        }
    }
 
TEST(TieIndices)
    {
    ITensor T(l1,l2,a1,s2,s1);
    T.Randomize();

    Index tied("tied",l2.m());
    ITSparse S(l2,l1,s1,tied,1);

    ITensor TT = S * T;

    CHECK_EQUAL(TT.r(),3);

    for(int j = 1; j <= 2; ++j)
    for(int k = 1; k <= 2; ++k)
        {
        CHECK_CLOSE(T(l1(j),l2(j),a1(1),s2(k),s1(j)),TT(tied(j),s2(k),a1(1)),1E-5);
        }
    }

TEST(Trace)
    {

    ITensor A(b2,a1,b3,b5,primed(b3));
    A.Randomize();
    Real f = -ran1();
    A *= f;

    ITensor At = ITSparse(b3,primed(b3),1) * A;

    Vector v(b3.m()); v = 1;
    ITensor AtV = ITSparse(b3,primed(b3),v) * A;

    for(int j2 = 1; j2 <= b2.m(); ++j2)
    for(int j5 = 1; j5 <= b5.m(); ++j5)
        {
        Real val = 0;
        for(int j3 = 1; j3 <= b3.m(); ++j3)
            {
            val += A(b2(j2),a1(1),b3(j3),b5(j5),primed(b3)(j3));
            }
        CHECK_CLOSE(val,At(b2(j2),a1(1),b5(j5)),1E-10);
        CHECK_CLOSE(val,AtV(b2(j2),a1(1),b5(j5)),1E-10);
        }
    }


BOOST_AUTO_TEST_SUITE_END()
