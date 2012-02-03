//#define THIS_IS_MAIN
#include "test.h"
#include "itensor.h"
#include <boost/test/unit_test.hpp>

struct ITensorDefaults
{
    Index s1,s2,s3,s4,
          s1P,s2P,s3P,s4P,
          l1,l2,l3,l4,l5,l6,l7,l8,
          a1,a2,a3,a4,
          b2,b3,b4,b5;

    ITensor A,B,X,Z;

    std::vector<Index> mixed_inds, reordered_mixed_inds;

    int mixed_inds_dim;

    ITensorDefaults() :
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

    ~ITensorDefaults() { }

};

BOOST_FIXTURE_TEST_SUITE(ITensorTest,ITensorDefaults)

TEST(Null)
{
    ITensor t1;

    CHECK(t1.isNull());

    ITensor t2(s1);

    CHECK(t2.isNotNull());
}

TEST(Constructors)
{
    ITensor t1(l1);

    CHECK_EQUAL(t1.r(),1);
    CHECK(t1.hasindex(l1));
    CHECK_CLOSE(t1.norm(),0,1E-10);

    ITensor t2(l1,l2);

    CHECK_EQUAL(t2.r(),2);
    CHECK(t2.hasindex(l1));
    CHECK(t2.hasindex(l2));
    CHECK_CLOSE(t2.norm(),0,1E-10);

    ITensor t3(l1,l2,l3);

    CHECK_EQUAL(t3.r(),3);
    CHECK(t3.hasindex(l1));
    CHECK(t3.hasindex(l2));
    CHECK(t3.hasindex(l3));
    CHECK_CLOSE(t3.norm(),0,1E-10);

    ITensor t4(a1,l1);

    CHECK_EQUAL(t4.r(),2);
    CHECK(t4.hasindex(a1));
    CHECK(t4.hasindex(l1));
    CHECK_CLOSE(t4.norm(),0,1E-10);

    ITensor t5(l1,a1,l2);

    CHECK_EQUAL(t5.r(),3);
    CHECK(t5.hasindex(a1));
    CHECK(t5.hasindex(l1));
    CHECK(t5.hasindex(l2));
    CHECK_CLOSE(t5.norm(),0,1E-10);

    ITensor t6(l1,a1,l2,a2);

    CHECK_EQUAL(t6.r(),4);
    CHECK(t6.hasindex(l1));
    CHECK(t6.hasindex(a1));
    CHECK(t6.hasindex(l2));
    CHECK(t6.hasindex(a2));
    CHECK_CLOSE(t6.norm(),0,1E-10);

    Real a = ran1();
    ITensor t7(l1,l2,a);

    CHECK_EQUAL(t7.r(),2);
    CHECK(t7.hasindex(l1));
    CHECK(t7.hasindex(l2));
    CHECK_CLOSE(t7(l1(1),l2(1)),a,1E-10);
    CHECK_CLOSE(t7(l1(1),l2(2)),0,1E-10);
    CHECK_CLOSE(t7(l1(2),l2(1)),0,1E-10);
    CHECK_CLOSE(t7(l1(2),l2(2)),a,1E-10);
    CHECK_CLOSE(t7.norm(),sqrt(min(l1.m(),l2.m()))*fabs(a),1E-10);

    Matrix M(l1.m(),l2.m()); 
    M(1,1) = ran1(); M(1,2) = ran1();
    M(2,1) = ran1(); M(2,2) = ran1();
    ITensor t8(l1,l2,M);

    CHECK_EQUAL(t8.r(),2);
    CHECK(t8.hasindex(l1));
    CHECK(t8.hasindex(l2));
    CHECK_CLOSE(t8(l1(1),l2(1)),M(1,1),1E-10);
    CHECK_CLOSE(t8(l1(1),l2(2)),M(1,2),1E-10);
    CHECK_CLOSE(t8(l1(2),l2(1)),M(2,1),1E-10);
    CHECK_CLOSE(t8(l1(2),l2(2)),M(2,2),1E-10);
    CHECK_CLOSE(t8.sumels(),M.TreatAsVector().sumels(),1E-10);
    CHECK_CLOSE(t8.norm(),Norm(M.TreatAsVector()),1E-10);

    Matrix W(a1.m(),l2.m()); 
    W(1,1) = ran1(); W(1,2) = ran1();
    ITensor w1(a1,l2,W);

    CHECK_EQUAL(w1.r(),2);
    CHECK(w1.hasindex(a1));
    CHECK(w1.hasindex(l2));
    CHECK_CLOSE(w1(l2(1)),W(1,1),1E-10);
    CHECK_CLOSE(w1(l2(2)),W(1,2),1E-10);
    CHECK_CLOSE(w1.sumels(),W.TreatAsVector().sumels(),1E-10);
    CHECK_CLOSE(w1.norm(),Norm(W.TreatAsVector()),1E-10);

    ITensor w2(l2,a1,W.t());

    CHECK_EQUAL(w2.r(),2);
    CHECK(w2.hasindex(a1));
    CHECK(w2.hasindex(l2));
    CHECK_CLOSE(w2(l2(1)),W(1,1),1E-10);
    CHECK_CLOSE(w2(l2(2)),W(1,2),1E-10);
    CHECK_CLOSE(w2.sumels(),W.TreatAsVector().sumels(),1E-10);
    CHECK_CLOSE(w2.norm(),Norm(W.TreatAsVector()),1E-10);

    Real b = ran1();
    ITensor t9(b);

    CHECK_CLOSE(t9.sumels(),b,1E-10);
    CHECK_CLOSE(t9.norm(),fabs(b),1E-10);

    Index linkind("linkind",10);
    Vector V(linkind.m()); V.Randomize();
    ITensor t10(linkind,V);

    CHECK_EQUAL(t10.r(),1);
    CHECK(t10.hasindex(linkind));
    CHECK_CLOSE(t10.sumels(),V.sumels(),1E-10);
    CHECK_CLOSE(t10.norm(),Norm(V),1E-10);
}

TEST(IndexValConstructors)
{
    ITensor t1(l1(2));

    CHECK_EQUAL(t1.r(),1);
    CHECK(t1.hasindex(l1));
    CHECK_CLOSE(t1(l1(1)),0,1E-10);
    CHECK_CLOSE(t1(l1(2)),1,1E-10);
    CHECK_CLOSE(t1.sumels(),1,1E-10);
    CHECK_CLOSE(t1.norm(),1,1E-10);

    ITensor t2(l1(2),l2(1));

    CHECK_EQUAL(t2.r(),2);
    CHECK(t2.hasindex(l1));
    CHECK(t2.hasindex(l2));
    CHECK_CLOSE(t2(l1(1),l2(1)),0,1E-10);
    CHECK_CLOSE(t2(l1(1),l2(2)),0,1E-10);
    CHECK_CLOSE(t2(l1(2),l2(1)),1,1E-10);
    CHECK_CLOSE(t2(l1(2),l2(2)),0,1E-10);
    CHECK_CLOSE(t2.sumels(),1,1E-10);
    CHECK_CLOSE(t2.norm(),1,1E-10);

    ITensor u2a(a1(1),l2(2));

    CHECK_EQUAL(u2a.r(),2);
    CHECK(u2a.hasindex(a1));
    CHECK(u2a.hasindex(l2));
    CHECK_CLOSE(u2a(l2(1)),0,1E-10);
    CHECK_CLOSE(u2a(l2(2)),1,1E-10);
    CHECK_CLOSE(u2a.sumels(),1,1E-10);
    CHECK_CLOSE(u2a.norm(),1,1E-10);

    ITensor u2b(l1(2),a2(1));

    CHECK_EQUAL(u2b.r(),2);
    CHECK(u2b.hasindex(l1));
    CHECK(u2b.hasindex(a2));
    CHECK_CLOSE(u2b(l1(1)),0,1E-10);
    CHECK_CLOSE(u2b(l1(2)),1,1E-10);
    CHECK_CLOSE(u2b.sumels(),1,1E-10);
    CHECK_CLOSE(u2b.norm(),1,1E-10);

    ITensor t3(l1(2),l3(1),l2(1));

    CHECK_EQUAL(t3.r(),3);
    CHECK(t3.hasindex(l1));
    CHECK(t3.hasindex(l2));
    CHECK(t3.hasindex(l3));
    CHECK_CLOSE(t3(l1(1),l3(1),l2(1)),0,1E-10);
    CHECK_CLOSE(t3(l1(2),l3(1),l2(1)),1,1E-10);
    CHECK_CLOSE(t3(l1(1),l3(2),l2(1)),0,1E-10);
    CHECK_CLOSE(t3(l1(2),l3(2),l2(2)),0,1E-10);
    CHECK_CLOSE(t3(l1(1),l3(1),l2(2)),0,1E-10);
    CHECK_CLOSE(t3(l1(2),l3(1),l2(2)),0,1E-10);
    CHECK_CLOSE(t3(l1(1),l3(2),l2(2)),0,1E-10);
    CHECK_CLOSE(t3(l1(2),l3(2),l2(2)),0,1E-10);
    CHECK_CLOSE(t3.sumels(),1,1E-10);
    CHECK_CLOSE(t3.norm(),1,1E-10);

    ITensor t4(a1(1),l3(2),l2(1));

    CHECK_EQUAL(t4.r(),3);
    CHECK(t4.hasindex(a1));
    CHECK(t4.hasindex(l2));
    CHECK(t4.hasindex(l3));
    CHECK_CLOSE(t4(l3(1),l2(1)),0,1E-10);
    CHECK_CLOSE(t4(l3(1),l2(2)),0,1E-10);
    CHECK_CLOSE(t4(l3(2),l2(1)),1,1E-10);
    CHECK_CLOSE(t4(l3(2),l2(2)),0,1E-10);
    CHECK_CLOSE(t4.sumels(),1,1E-10);
    CHECK_CLOSE(t4.norm(),1,1E-10);
}

TEST(MultiIndexConstructors)
{
    std::vector<Index> indices(4);
    indices[0] = a2;
    indices[1] = l3;
    indices[2] = l1;
    indices[3] = a4;

    ITensor t1(indices);

    CHECK_EQUAL(t1.r(),4);
    CHECK(t1.hasindex(a2));
    CHECK(t1.hasindex(l3));
    CHECK(t1.hasindex(l1));
    CHECK(t1.hasindex(a4));
    CHECK_CLOSE(t1.norm(),0,1E-10);

    Vector V(l1.m()*l3.m());
    V.Randomize();

    ITensor t2(indices,V);

    CHECK_EQUAL(t2.r(),4);
    CHECK(t2.hasindex(a2));
    CHECK(t2.hasindex(l3));
    CHECK(t2.hasindex(l1));
    CHECK(t2.hasindex(a4));
    CHECK_CLOSE(t2.norm(),Norm(V),1E-10);
    CHECK_CLOSE(t2.sumels(),V.sumels(),1E-10);
}

TEST(ITensorConstructors)
{
    Index clink("clink",4);
    std::vector<Index> indices1(3);
    indices1.at(0) = l1;
    indices1.at(1) = l2;
    indices1.at(2) = clink;

    Vector V(l1.m()*l2.m()*clink.m());
    V.Randomize();

    ITensor t1(indices1,V);

    Real f = ran1();

    ITensor t2(t1);
    t2 *= f;

    std::vector<Index> indices3(4);
    indices3.at(0) = l1;
    indices3.at(1) = l2;
    indices3.at(2) = l3;
    indices3.at(3) = l4;

    ITensor t3(indices3,t2);

    CHECK_EQUAL(4,t3.r());

    for(int i = 1; i <= l1.m(); ++i)
    for(int j = 1; j <= l2.m(); ++j)
    {
    CHECK_CLOSE(t1(l1(i),l2(j),clink(1))*f,t3(l1(i),l2(j),l3(1),l4(1)),1E-10);
    CHECK_CLOSE(t1(l1(i),l2(j),clink(2))*f,t3(l1(i),l2(j),l3(2),l4(1)),1E-10);
    CHECK_CLOSE(t1(l1(i),l2(j),clink(3))*f,t3(l1(i),l2(j),l3(1),l4(2)),1E-10);
    CHECK_CLOSE(t1(l1(i),l2(j),clink(4))*f,t3(l1(i),l2(j),l3(2),l4(2)),1E-10);
    }

    Permutation P;
    P.from_to(2,4);
    P.from_to(4,2);
    CHECK(P.check(4));

    std::vector<Index> indices5(4);
    indices5.at(0) = l1;
    indices5.at(1) = l4;
    indices5.at(2) = l3;
    indices5.at(3) = l2;

    ITensor t4(t3);
    Real f2 = ran1();
    t4 /= f2;
    ITensor t5(indices5,t4,P);

    CHECK_EQUAL(4,t5.r());

    for(int i = 1; i <= l1.m(); ++i)
    for(int j = 1; j <= l2.m(); ++j)
    for(int k = 1; k <= l3.m(); ++k)
    for(int l = 1; l <= l4.m(); ++l)
    {
    CHECK_CLOSE(t3(l1(i),l2(j),l3(k),l4(l))/f2,t5(l1(i),l2(j),l3(k),l4(l)),1E-10);
    }

}

TEST(Copy)
{
    std::vector<Index> indices(4);
    indices[0] = a2;
    indices[1] = l3;
    indices[2] = l1;
    indices[3] = a4;

    Vector V(l1.m()*l3.m());
    V.Randomize();

    ITensor t1(indices,V);

    CHECK_EQUAL(t1.r(),4);
    CHECK(t1.hasindex(a2));
    CHECK(t1.hasindex(l3));
    CHECK(t1.hasindex(l1));
    CHECK(t1.hasindex(a4));
    CHECK_CLOSE(t1.norm(),Norm(V),1E-10);
    CHECK_CLOSE(t1.sumels(),V.sumels(),1E-10);

    //Use copy constructor
    ITensor t2(t1);
    t1 = ITensor(); //destroy t1

    CHECK_EQUAL(t2.r(),4);
    CHECK(t2.hasindex(a2));
    CHECK(t2.hasindex(l3));
    CHECK(t2.hasindex(l1));
    CHECK(t2.hasindex(a4));
    CHECK_CLOSE(t2.norm(),Norm(V),1E-10);
    CHECK_CLOSE(t2.sumels(),V.sumels(),1E-10);

    //Use operator=
    ITensor t3 = t2;
    t2 = ITensor(); //destroy t2

    CHECK_EQUAL(t3.r(),4);
    CHECK(t3.hasindex(a2));
    CHECK(t3.hasindex(l3));
    CHECK(t3.hasindex(l1));
    CHECK(t3.hasindex(a4));
    CHECK_CLOSE(t3.norm(),Norm(V),1E-10);
    CHECK_CLOSE(t3.sumels(),V.sumels(),1E-10);
}

TEST(ScalarMultiply)
{
    A *= -1;
    CHECK_EQUAL(A(s1(1),s2(1)),-11);
    CHECK_EQUAL(A(s1(1),s2(2)),-12);
    CHECK_EQUAL(A(s1(2),s2(1)),-21);
    CHECK_EQUAL(A(s1(2),s2(2)),-22);

    Real f = ran1();
    A *= -f;
    CHECK_CLOSE(A(s1(1),s2(1)),11*f,1E-10);
    CHECK_CLOSE(A(s1(1),s2(2)),12*f,1E-10);
    CHECK_CLOSE(A(s1(2),s2(1)),21*f,1E-10);
    CHECK_CLOSE(A(s1(2),s2(2)),22*f,1E-10);

    B /= f;
    CHECK_CLOSE(B(s1(1),s2(1)),110/f,1E-10);
    CHECK_CLOSE(B(s1(1),s2(2)),120/f,1E-10);
    CHECK_CLOSE(B(s1(2),s2(1)),210/f,1E-10);
    CHECK_CLOSE(B(s1(2),s2(2)),220/f,1E-10);
}

TEST(assignToVec)
{
    Vector V(l1.m()*l2.m()*l3.m());
    V.Randomize();
    Real f = -ran1();

    std::vector<Index> indices; indices.reserve(3);
    indices.push_back(l1);
    indices.push_back(l2);
    indices.push_back(l3);

    ITensor T(indices,V);

    T *= f;

    Vector U(T.vecSize()); T.assignToVec(U);

    CHECK_EQUAL(U.Length(),V.Length());

    for(int j = 1; j < V.Length(); ++j)
    { CHECK_CLOSE(U(j),V(j)*f,1E-10); }

}

TEST(MapElems)
    {
    // class Functor and the function Func
    // are defined in test.h

    ITensor A1(A);
    Functor f;
    A1.mapElems(f);
    for(int n1 = 1; n1 <= s1.m(); ++n1)
    for(int n2 = 1; n2 <= s2.m(); ++n2)
        {
        CHECK_CLOSE( f( A(s1(n1),s2(n2)) ), A1(s1(n1),s2(n2)) ,1E-10);
        }

    ITensor A2(A);
    Real (*pFunc)(Real) = Func;
    A2.mapElems(*pFunc);
    for(int n1 = 1; n1 <= s1.m(); ++n1)
    for(int n2 = 1; n2 <= s2.m(); ++n2)
        {
        CHECK_CLOSE( Func( A(s1(n1),s2(n2)) ), A2(s1(n1),s2(n2)) ,1E-10);
        }
    }

TEST(reshape)
{
    Permutation P;
    P.from_to(1,2);
    P.from_to(2,1);

    Real f = -5;
    A *= f;

    A.reshape(P);

    CHECK_CLOSE(A(s1(1),s2(1)),11*f,1E-10);
    CHECK_CLOSE(A(s1(1),s2(2)),21*f,1E-10);
    CHECK_CLOSE(A(s1(2),s2(1)),12*f,1E-10);
    CHECK_CLOSE(A(s1(2),s2(2)),22*f,1E-10);

}

TEST(findindex)
{
    ITensor T(mixed_inds);

    boost::array<int,6> arb_order = {{ 3, 4, 1, 0, 2, 5 }};

    Foreach(int i, arb_order)
    {
        int j = T.findindex(mixed_inds.at(i));
        CHECK_EQUAL(T.index(j),mixed_inds.at(i));
    }

}

TEST(SumDifference)
{
    Vector V(mixed_inds_dim),W(mixed_inds_dim);
    V.Randomize();
    W.Randomize();

    ITensor v(mixed_inds,V), w(mixed_inds,W);

    Real f1 = -ran1(), f2 = 0.1*f1;

    ITensor r = f1*v + w/f2; 
    Vector R(r.vecSize()); r.assignToVec(R);
    for(int j = 1; j < R.Length(); ++j)
    { CHECK_CLOSE(R(j),f1*V(j)+W(j)/f2,1E-10); }

    ITensor d(v); d -= w;
    Vector D(d.vecSize()); d.assignToVec(D);
    for(int j = 1; j < D.Length(); ++j)
    { CHECK_CLOSE(D(j),V(j)-W(j),1E-10); }

    Vector YY(mixed_inds_dim),ZZ(mixed_inds_dim);
    YY.Randomize();
    ZZ.Randomize();

    f1 = 1; f2 = 1;
    ITensor yy(mixed_inds,YY), zz(reordered_mixed_inds,ZZ);
    r = f1*yy + f2*zz;
    for(int j1 = 1; j1 <= 2; ++j1)
    for(int j2 = 1; j2 <= 2; ++j2)
    for(int k3 = 1; k3 <= 3; ++k3)
    for(int j4 = 1; j4 <= 2; ++j4)
    { 
        CHECK_CLOSE(r(l1(j1),l2(j2),b3(k3),l4(j4)),
                    f1*yy(l1(j1),l2(j2),b3(k3),l4(j4))
                   +f2*zz(l1(j1),l2(j2),b3(k3),l4(j4))
                   ,1E-10); 
    }

}

TEST(ContractingProduct)
{

    //Check for rank 0 ITensors
    {
    Real f = ran1();
    ITensor rZ(f), T(b2,a1,b4);
    T.Randomize();

    ITensor res = rZ * T;

    CHECK_EQUAL(rZ.r(),0);
    CHECK_EQUAL(res.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
        {
        Real val = f * T(b2(j2),a1(1),b4(j4));
        CHECK_CLOSE(res(b2(j2),a1(1),b4(j4)),val,1E-10);
        }
    }

    //More general case
    ITensor L(b2,a1,b3,b4), R(a1,b3,a2,b5,b4);

    L.Randomize(); R.Randomize();

    Real fL = ran1(), fR = ran1();
    ITensor Lf = L * fL;
    ITensor Rf = R * fR;

    ITensor res1 = Lf*Rf;

    CHECK(res1.hasindex(b2));
    CHECK(res1.hasindex(a2));
    CHECK(res1.hasindex(b5));
    CHECK(!res1.hasindex(a1));
    CHECK(!res1.hasindex(b3));
    CHECK(!res1.hasindex(b4));

    CHECK_EQUAL(res1.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j5 = 1; j5 <= 5; ++j5)
    {
        Real val = 0;
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int j4 = 1; j4 <= 4; ++j4)
        {
            val += L(b2(j2),a1(1),b3(j3),b4(j4))*fL * R(a1(1),b3(j3),a2(1),b5(j5),b4(j4))*fR;
        }
        CHECK_CLOSE(res1(b2(j2),a2(1),b5(j5)),val,1E-10);
    }

    ITensor res2 = R*L;

    CHECK(res2.hasindex(b2));
    CHECK(res2.hasindex(a2));
    CHECK(res2.hasindex(b5));
    CHECK(!res2.hasindex(a1));
    CHECK(!res2.hasindex(b3));
    CHECK(!res2.hasindex(b4));

    CHECK_EQUAL(res2.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j5 = 1; j5 <= 5; ++j5)
    {
        Real val = 0;
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int j4 = 1; j4 <= 4; ++j4)
        {
            val += L(b2(j2),a1(1),b3(j3),b4(j4)) * R(a1(1),b3(j3),a2(1),b5(j5),b4(j4));
        }
        CHECK_CLOSE(res2(b2(j2),a2(1),b5(j5)),val,1E-10);
    }

    ITensor Q(a1,b4,a2,b2), P(a2,a3,a1);

    Q.Randomize(); P.Randomize();

    Real fQ = ran1(), fP = ran1();
    ITensor Qf = Q * fQ;
    ITensor Pf = P * fP;

    ITensor res3 = Qf*Pf;

    CHECK(res3.hasindex(b4));
    CHECK(res3.hasindex(b2));
    CHECK(res3.hasindex(a3));
    CHECK(!res3.hasindex(a1));
    CHECK(!res3.hasindex(a2));

    CHECK_EQUAL(res3.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
    {
        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
        CHECK_CLOSE(res3(b4(j4),b2(j2)),val,1E-10);
    }

    ITensor res4 = Pf*Qf;

    CHECK(res4.hasindex(b4));
    CHECK(res4.hasindex(b2));
    CHECK(res4.hasindex(a3));
    CHECK(!res4.hasindex(a1));
    CHECK(!res4.hasindex(a2));

    CHECK_EQUAL(res4.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
    {
        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
        CHECK_CLOSE(res4(b4(j4),b2(j2)),val,1E-10);
    }


    ITensor psi(a1,a2,a3), mpoh(l2,a1,a1.primed(),a2,a2.primed());
    psi.Randomize(); mpoh.Randomize();

    ITensor Hpsi = mpoh * psi;

    CHECK_EQUAL(Hpsi.r(),4);
    CHECK(Hpsi.hasindex(l2));
    CHECK(Hpsi.hasindex(a1.primed()));
    CHECK(Hpsi.hasindex(a2.primed()));
    CHECK(Hpsi.hasindex(a3));
    CHECK(!Hpsi.hasindex(a1));
    CHECK(!Hpsi.hasindex(a2));
}

TEST(NonContractingProduct)
{
    ITensor L(b2,a1,b3,b4), R(a1,b3,a2,b5,b4);

    L.Randomize(); R.Randomize();

    Real fL = ran1(), fR = ran1();
    ITensor Lf = L * fL;
    ITensor Rf = R * fR;

    ITensor res1 = Lf / Rf;

    CHECK(res1.hasindex(b2));
    CHECK(res1.hasindex(a2));
    CHECK(res1.hasindex(b5));
    CHECK(res1.hasindex(a1));
    CHECK(res1.hasindex(b3));
    CHECK(res1.hasindex(b4));

    CHECK_EQUAL(res1.r(),6);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j3 = 1; j3 <= 3; ++j3)
    for(int j4 = 1; j4 <= 4; ++j4)
    for(int j5 = 1; j5 <= 5; ++j5)
    {
        Real val = L(b2(j2),b3(j3),b4(j4))*fL * R(b3(j3),b5(j5),b4(j4))*fR;
        CHECK_CLOSE(res1(b2(j2),b3(j3),b4(j4),b5(j5)),val,1E-10);
    }

    ITensor res2 = R/L;

    CHECK(res2.hasindex(b2));
    CHECK(res2.hasindex(a2));
    CHECK(res2.hasindex(b5));
    CHECK(res2.hasindex(a1));
    CHECK(res2.hasindex(b3));
    CHECK(res2.hasindex(b4));

    CHECK_EQUAL(res2.r(),6);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j3 = 1; j3 <= 3; ++j3)
    for(int j4 = 1; j4 <= 4; ++j4)
    for(int j5 = 1; j5 <= 5; ++j5)
    {
        Real val = L(b2(j2),a1(1),b3(j3),b4(j4)) * R(a1(1),b3(j3),a2(1),b5(j5),b4(j4));
        CHECK_CLOSE(res2(b2(j2),b3(j3),b4(j4),b5(j5)),val,1E-10);
    }

    ITensor Q(a1,b4,a2,b2), P(a2,a3,a1);

    Q.Randomize(); P.Randomize();

    Real fQ = ran1(), fP = ran1();
    ITensor Qf = Q * fQ;
    ITensor Pf = P * fP;

    ITensor res3 = Qf/Pf;

    CHECK(res3.hasindex(b4));
    CHECK(res3.hasindex(b2));
    CHECK(res3.hasindex(a3));
    CHECK(res3.hasindex(a1));
    CHECK(res3.hasindex(a2));

    CHECK_EQUAL(res3.r(),5);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
    {
        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
        CHECK_CLOSE(res3(b4(j4),b2(j2)),val,1E-10);
    }

    ITensor res4 = Pf/Qf;

    CHECK(res4.hasindex(b4));
    CHECK(res4.hasindex(b2));
    CHECK(res4.hasindex(a3));
    CHECK(res4.hasindex(a1));
    CHECK(res4.hasindex(a2));

    CHECK_EQUAL(res4.r(),5);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
    {
        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
        CHECK_CLOSE(res4(b4(j4),b2(j2)),val,1E-10);
    }


    ITensor psi(a1,a2,a3), mpoh(l2,a1,a1.primed(),a2,a2.primed());
    psi.Randomize(); mpoh.Randomize();

    ITensor Hpsi = mpoh / psi;

    CHECK_EQUAL(Hpsi.r(),6);
    CHECK(Hpsi.hasindex(l2));
    CHECK(Hpsi.hasindex(a1));
    CHECK(Hpsi.hasindex(a2));
    CHECK(Hpsi.hasindex(a1.primed()));
    CHECK(Hpsi.hasindex(a2.primed()));
    CHECK(Hpsi.hasindex(a3));

    for(int j2 = 1; j2 <= 2; ++j2)
    { CHECK_CLOSE(Hpsi(l2(j2)),psi()*mpoh(l2(j2)),1E-10); }
}

TEST(TieIndices)
    {

    Index t("tied",2);

    ITensor dX(X);
    dX.tieIndices(s1,s2,t);

    CHECK_CLOSE(dX.norm(),0,1E-5);
    CHECK_EQUAL(dX.r(),1);
    CHECK(dX.hasindex(t));

    ITensor dZ(Z);
    dZ.tieIndices(s1,s2,t);
    CHECK_CLOSE(dZ(t(1)),+1,1E-5);
    CHECK_CLOSE(dZ(t(2)),-1,1E-5);

    {
    ITensor T(l1,l2,a1,s2,s1);
    T.Randomize();

    ITensor TT(T);
    TT.tieIndices(l2,l1,s1,l2);

    CHECK_EQUAL(TT.r(),3);

    for(int j = 1; j <= 2; ++j)
    for(int k = 1; k <= 2; ++k)
        {
        CHECK_CLOSE(T(l1(j),l2(j),a1(1),s2(k),s1(j)),TT(l2(j),s2(k),a1(1)),1E-5);
        }
    }

    //Try tying m==1 inds
    {
    ITensor T(l1,a2,a1,s2,a3);
    T.Randomize();

    ITensor TT(T);
    TT.tieIndices(a1,a3,a2,a1);

    CHECK_EQUAL(TT.r(),3);

    for(int j = 1; j <= 2; ++j)
    for(int k = 1; k <= 2; ++k)
        {
        CHECK_CLOSE(T(l1(j),s2(k)),TT(l1(j),s2(k)),1E-5);
        }
    }



    } //TieIndices

TEST(fromMatrix11)
{
    Matrix M22(s1.m(),s2.m());

    M22(1,1) = -0.3; M22(1,2) = 110;
    M22(2,1) = -1.7; M22(1,2) = 5;

    ITensor T(s1,s2);
    //T should be overwritten so check
    //that scalar mult has no effect
    T *= -5; 

    T.fromMatrix11(s1,s2,M22);

    CHECK_CLOSE(T(s1(1),s2(1)),M22(1,1),1E-10);
    CHECK_CLOSE(T(s1(1),s2(2)),M22(1,2),1E-10);
    CHECK_CLOSE(T(s1(2),s2(1)),M22(2,1),1E-10);
    CHECK_CLOSE(T(s1(2),s2(2)),M22(2,2),1E-10);

    ITensor U(T);

    U.fromMatrix11(s2,s1,M22);

    CHECK_CLOSE(T(s1(1),s2(1)),M22(1,1),1E-10);
    CHECK_CLOSE(T(s1(1),s2(2)),M22(1,2),1E-10);
    CHECK_CLOSE(T(s1(2),s2(1)),M22(2,1),1E-10);
    CHECK_CLOSE(T(s1(2),s2(2)),M22(2,2),1E-10);

    CHECK_CLOSE(U(s2(1),s1(1)),M22(1,1),1E-10);
    CHECK_CLOSE(U(s2(1),s1(2)),M22(1,2),1E-10);
    CHECK_CLOSE(U(s2(2),s1(1)),M22(2,1),1E-10);
    CHECK_CLOSE(U(s2(2),s1(2)),M22(2,2),1E-10);

    Matrix M12(a1.m(),s2.m());
    M12(1,1) = 37; M12(1,2) = -2;

    ITensor P(a1,s2);
    P *= -4;

    P.fromMatrix11(a1,s2,M12);

    CHECK_CLOSE(P(a1(1),s2(1)),M12(1,1),1E-10);
    CHECK_CLOSE(P(a1(1),s2(2)),M12(1,2),1E-10);

    P.fromMatrix11(s2,a1,M12.t());

    CHECK_CLOSE(P(s2(1),a1(1)),M12(1,1),1E-10);
    CHECK_CLOSE(P(s2(2),a1(1)),M12(1,2),1E-10);
}

TEST(toMatrix11)
    {
    Matrix M(s1.m(),s2.m());    

    Real f = -ran1();

    A *= f;

    A.toMatrix11(s2,s1,M);

    CHECK_CLOSE(M(1,1),11*f,1E-10);
    CHECK_CLOSE(M(2,1),12*f,1E-10);
    CHECK_CLOSE(M(1,2),21*f,1E-10);
    CHECK_CLOSE(M(2,2),22*f,1E-10);

    A.toMatrix11(s1,s2,M);

    CHECK_CLOSE(M(1,1),11*f,1E-10);
    CHECK_CLOSE(M(1,2),12*f,1E-10);
    CHECK_CLOSE(M(2,1),21*f,1E-10);
    CHECK_CLOSE(M(2,2),22*f,1E-10);

    A.toMatrix11NoScale(s2,s1,M);

    CHECK_CLOSE(M(1,1),11,1E-10);
    CHECK_CLOSE(M(2,1),12,1E-10);
    CHECK_CLOSE(M(1,2),21,1E-10);
    CHECK_CLOSE(M(2,2),22,1E-10);

    Vector V(4);
    V(1) = 3.14; V(2) = 2.718; V(3) = -1; V(4) = 0;
    Index link("link",4);

    ITensor T(link,a1);
    T.assignFromVec(V);

    Matrix M41(4,1), M14(1,4);
    
    T.toMatrix11(link,a1,M41);

    CHECK_CLOSE(M41(1,1),V(1),1E-10);
    CHECK_CLOSE(M41(2,1),V(2),1E-10);
    CHECK_CLOSE(M41(3,1),V(3),1E-10);
    CHECK_CLOSE(M41(4,1),V(4),1E-10);
     
    T.toMatrix11(a1,link,M14);

    CHECK_CLOSE(M14(1,1),V(1),1E-10);
    CHECK_CLOSE(M14(1,2),V(2),1E-10);
    CHECK_CLOSE(M14(1,3),V(3),1E-10);
    CHECK_CLOSE(M14(1,4),V(4),1E-10);

    }

TEST(SymmetricDiag11)
    {
    ITensor T(s1,s1.primed());
    commaInit(T) << 1, 2,
                    2, 1;
    T *= -2;
    Index mid;
    ITensor D,U;
    T.symmetricDiag11(s1,D,U,mid);
    ITensor UD(U);
    UD.primeind(s1);
    UD /= D;
    ITensor diff(UD*U - T);
    CHECK(diff.norm() < 1E-10);

    //Construct a random, symmetric ITensor
    const int qs = 50;
    Index q("q",qs);
    Matrix QQ(qs,qs);
    QQ.Randomize();
    QQ += QQ.t();
    ITensor Q(q,primed(q),QQ);
    Q *= -2;

    //Diagonalize and check the factorization
    Q.symmetricDiag11(q,D,U,mid);
    UD =U;
    UD.primeind(q);
    UD /= D;
    diff = UD*U - Q;
    CHECK(diff.norm() < 1E-10);
    }

TEST(CommaAssignment)
    {
    ITensor ZZ(s1,s2);
    commaInit(ZZ) << 1, 0, 
                     0, -1;
    CHECK_EQUAL(ZZ(s1(1),s2(1)),1);
    CHECK_EQUAL(ZZ(s1(2),s2(1)),0);
    CHECK_EQUAL(ZZ(s1(1),s2(2)),0);
    CHECK_EQUAL(ZZ(s1(2),s2(2)),-1);

    ITensor XX(s1,s2);
    XX(s1(2),s2(1)) = 5;
    XX *= 3;
    commaInit(XX) << 0, 1, 
                     1, 0;
    CHECK_EQUAL(XX(s1(1),s2(1)),0);
    CHECK_EQUAL(XX(s1(2),s2(1)),1);
    CHECK_EQUAL(XX(s1(1),s2(2)),1);
    CHECK_EQUAL(XX(s1(2),s2(2)),0);

    ITensor AA(s1,s2);
    AA.Randomize();
    AA *= -ran1();
    commaInit(AA) << 11, 21, 
                     12, 22;
    CHECK_EQUAL(AA(s1(1),s2(1)),11);
    CHECK_EQUAL(AA(s1(2),s2(1)),21);
    CHECK_EQUAL(AA(s1(1),s2(2)),12);
    CHECK_EQUAL(AA(s1(2),s2(2)),22);
    }

TEST(Website)
    {

    Index a("a",2), b("b",2), c("c",2);
    ITensor Z(a,b), X(b,c);
    commaInit(Z) << 1, 0, 
                    0, -1;
    commaInit(X) << 0, 1, 
                    1, 0;
    ITensor R = Z * X;

    CHECK_CLOSE(R(a(1),c(1)),0,1E-10);
    CHECK_CLOSE(R(a(1),c(2)),+1,1E-10);
    CHECK_CLOSE(R(a(2),c(1)),-1,1E-10);
    CHECK_CLOSE(R(a(2),c(2)),0,1E-10);

    }

BOOST_AUTO_TEST_SUITE_END()
