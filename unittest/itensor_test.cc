#include "test.h"
#include "itensor.h"

using namespace std;
using namespace itensor;

//double static
//Func(double x)
//    {
//    return x*x;
//    }

class Functor
    {
    public:

    double
    operator()(double x) const
        {
        return x*x;
        }

    int
    operator()(int x) const
        {
        return x*x;
        }

    };

TEST_CASE("ITensor")
{

    Index s1("s1",2,Site);
    Index s2("s2",2,Site);
    Index s3("s3",2,Site);
    Index s4("s4",2,Site);
    Index s1P(prime(s1));
    Index s2P(prime(s2));
    Index s3P(prime(s3));
    Index s4P(prime(s4));
    Index l1("l1",2);
    Index l2("l2",2);
    Index l3("l3",2);
    Index l4("l4",2);
    Index l5("l5",2);
    Index l6("l6",2);
    Index l7("l7",2);
    Index l8("l8",2);
    Index a1("a1");
    Index a2("a2");
    Index a3("a3");
    Index a4("a4");
    Index b2("b2",2);
    Index b3("b3",3);
    Index b4("b4",4);
    Index b5("b5",5);

    IndexSet<Index> mixed_inds(a2,b3,l1,l2,a4,l4);
    //const int mixed_inds_dim = mixed_inds.dim();

    ITensor A,
            B,
            X,
            Z;

     //   {
     //   Matrix M(s1.m(),s2.m());
     //   M(1,1) = 11; M(1,2) = 12;
     //   M(2,1) = 21; M(2,2) = 22;
     //   A = ITensor(s1,s2,M);
     //   }

     //   {
     //   Matrix M(s1.m(),s2.m());
     //   M(1,1) = 110; M(1,2) = 120;
     //   M(2,1) = 210; M(2,2) = 220;
     //   B = ITensor(s1,s2,M);
     //   }

     //   {
     //   Matrix M(s1.m(),s2.m());
     //   M(1,1) = 0; M(1,2) = 1;
     //   M(2,1) = 1; M(2,2) = 0;
     //   X = ITensor(s1,s2,M);
     //   }
     //   
     //   {
     //   Matrix M(s1.m(),s2.m());
     //   M(1,1) = 1; M(1,2) =  0;
     //   M(2,1) = 0; M(2,2) = -1;
     //   Z = ITensor(s1,s2,M);
     //   }


SECTION("Boolean")
    {
    ITensor t1;

    CHECK(!t1);

    ITensor t2(s1);

    CHECK(t2);
    }

SECTION("Constructors")
    {
    ITensor t1(l1);

    CHECK_EQUAL(t1.r(),1);
    CHECK(hasindex(t1,l1));
    CHECK_CLOSE(norm(t1),0,1E-10);

    ITensor t2(l1,l2);

    CHECK_EQUAL(t2.r(),2);
    CHECK(hasindex(t2,l1));
    CHECK(hasindex(t2,l2));
    CHECK_CLOSE(norm(t2),0,1E-10);

    ITensor t3(l1,l2,l3);

    CHECK_EQUAL(t3.r(),3);
    CHECK(hasindex(t3,l1));
    CHECK(hasindex(t3,l2));
    CHECK(hasindex(t3,l3));
    CHECK_CLOSE(norm(t3),0,1E-10);

    ITensor t4(a1,l1);

    CHECK_EQUAL(t4.r(),2);
    CHECK(hasindex(t4,a1));
    CHECK(hasindex(t4,l1));
    CHECK_CLOSE(norm(t4),0,1E-10);

    ITensor t5(l1,a1,l2);

    CHECK_EQUAL(t5.r(),3);
    CHECK(hasindex(t5,a1));
    CHECK(hasindex(t5,l1));
    CHECK(hasindex(t5,l2));
    CHECK_CLOSE(norm(t5),0,1E-10);

    ITensor t6(l1,a1,l2,a2);

    CHECK_EQUAL(t6.r(),4);
    CHECK(hasindex(t6,l1));
    CHECK(hasindex(t6,a1));
    CHECK(hasindex(t6,l2));
    CHECK(hasindex(t6,a2));
    CHECK_CLOSE(norm(t6),0,1E-10);

    ITensor t7(l1,l2);
    Real a = -0.83;
    t7.fill(a);

    CHECK_EQUAL(t7.r(),2);
    CHECK(hasindex(t7,l1));
    CHECK(hasindex(t7,l2));
    CHECK_CLOSE(t7.real(l1(1),l2(1)),a,1E-5);
    CHECK_CLOSE(t7.real(l1(1),l2(2)),a,1E-5);
    CHECK_CLOSE(t7.real(l1(2),l2(1)),a,1E-5);
    CHECK_CLOSE(t7.real(l1(2),l2(2)),a,1E-5);
    PrintData(t7);
//
//    Matrix M(l1.m(),b3.m()); 
//    M(1,1) = 11; M(1,2) = 12; M(1,3) = 13;
//    M(2,1) = 21; M(2,2) = 22; M(2,3) = 23;
//    ITensor t8(l1,b3,M);
//
//    CHECK_EQUAL(t8.r(),2);
//    CHECK(hasindex(t8,l1));
//    CHECK(hasindex(t8,b3));
//    CHECK_CLOSE(t8(l1(1),b3(1)),11,1E-10);
//    CHECK_CLOSE(t8(l1(1),b3(2)),12,1E-10);
//    CHECK_CLOSE(t8(l1(1),b3(3)),13,1E-10);
//    CHECK_CLOSE(t8(l1(2),b3(1)),21,1E-10);
//    CHECK_CLOSE(t8(l1(2),b3(2)),22,1E-10);
//    CHECK_CLOSE(t8(l1(2),b3(3)),23,1E-10);
//    CHECK_CLOSE(sumels(t8),M.TreatAsVector().sumels(),1E-10);
//    CHECK_CLOSE(t8.norm(),Norm(M.TreatAsVector()),1E-10);
//
//    ITensor t85(b3,l1,M.t());
//
//    CHECK_EQUAL(t85.r(),2);
//    CHECK(hasindex(t85,l1));
//    CHECK(hasindex(t85,b3));
//    CHECK_CLOSE(t85(l1(1),b3(1)),11,1E-10);
//    CHECK_CLOSE(t85(l1(1),b3(2)),12,1E-10);
//    CHECK_CLOSE(t85(l1(1),b3(3)),13,1E-10);
//    CHECK_CLOSE(t85(l1(2),b3(1)),21,1E-10);
//    CHECK_CLOSE(t85(l1(2),b3(2)),22,1E-10);
//    CHECK_CLOSE(t85(l1(2),b3(3)),23,1E-10);
//    CHECK_CLOSE(sumels(t85),M.TreatAsVector().sumels(),1E-10);
//    CHECK_CLOSE(t85.norm(),Norm(M.TreatAsVector()),1E-10);
//
//    Matrix W(a1.m(),l2.m()); 
//    W(1,1) = Global::random(); W(1,2) = Global::random();
//    ITensor w1(a1,l2,W);
//
//    CHECK_EQUAL(w1.r(),2);
//    CHECK(hasindex(w1,a1));
//    CHECK(hasindex(w1,l2));
//    CHECK_CLOSE(w1(l2(1)),W(1,1),1E-10);
//    CHECK_CLOSE(w1(l2(2)),W(1,2),1E-10);
//    CHECK_CLOSE(sumels(w1),W.TreatAsVector().sumels(),1E-10);
//    CHECK_CLOSE(w1.norm(),Norm(W.TreatAsVector()),1E-10);
//
//    ITensor w2(l2,a1,W.t());
//
//    CHECK_EQUAL(w2.r(),2);
//    CHECK(hasindex(w2,a1));
//    CHECK(hasindex(w2,l2));
//    CHECK_CLOSE(w2(l2(1)),W(1,1),1E-10);
//    CHECK_CLOSE(w2(l2(2)),W(1,2),1E-10);
//    CHECK_CLOSE(sumels(w2),W.TreatAsVector().sumels(),1E-10);
//    CHECK_CLOSE(w2.norm(),Norm(W.TreatAsVector()),1E-10);
//
//    Real b = Global::random();
//    ITensor t9(b);
//
//    CHECK_CLOSE(sumels(t9),b,1E-10);
//    CHECK_CLOSE(t9.norm(),fabs(b),1E-10);
//
//    Index linkind("linkind",10);
//    Vector V(linkind.m()); V.Randomize();
//    ITensor t10(linkind,V);
//
//    CHECK_EQUAL(t10.r(),1);
//    CHECK(hasindex(t10,linkind));
//    CHECK_CLOSE(sumels(t10),V.sumels(),1E-10);
//    CHECK_CLOSE(t10.norm(),Norm(V),1E-10);
    }
//
//SECTION("IndexValConstructors")
//    {
//    ITensor t1(l1(2));
//
//    CHECK_EQUAL(t1.r(),1);
//    CHECK(hasindex(t1,l1));
//    CHECK_CLOSE(t1(l1(1)),0,1E-10);
//    CHECK_CLOSE(t1(l1(2)),1,1E-10);
//    CHECK_CLOSE(sumels(t1),1,1E-10);
//    CHECK_CLOSE(t1.norm(),1,1E-10);
//
//    ITensor t2(l1(2),l2(1));
//
//    CHECK_EQUAL(t2.r(),2);
//    CHECK(hasindex(t2,l1));
//    CHECK(hasindex(t2,l2));
//    CHECK_CLOSE(t2(l1(1),l2(1)),0,1E-10);
//    CHECK_CLOSE(t2(l1(1),l2(2)),0,1E-10);
//    CHECK_CLOSE(t2(l1(2),l2(1)),1,1E-10);
//    CHECK_CLOSE(t2(l1(2),l2(2)),0,1E-10);
//    CHECK_CLOSE(sumels(t2),1,1E-10);
//    CHECK_CLOSE(t2.norm(),1,1E-10);
//
//    ITensor u2a(a1(1),l2(2));
//
//    CHECK_EQUAL(u2a.r(),2);
//    CHECK(hasindex(u2a,a1));
//    CHECK(hasindex(u2a,l2));
//    CHECK_CLOSE(u2a(l2(1)),0,1E-10);
//    CHECK_CLOSE(u2a(l2(2)),1,1E-10);
//    CHECK_CLOSE(sumels(u2a),1,1E-10);
//    CHECK_CLOSE(u2a.norm(),1,1E-10);
//
//    ITensor u2b(l1(2),a2(1));
//
//    CHECK_EQUAL(u2b.r(),2);
//    CHECK(hasindex(u2b,l1));
//    CHECK(hasindex(u2b,a2));
//    CHECK_CLOSE(u2b(l1(1)),0,1E-10);
//    CHECK_CLOSE(u2b(l1(2)),1,1E-10);
//    CHECK_CLOSE(sumels(u2b),1,1E-10);
//    CHECK_CLOSE(u2b.norm(),1,1E-10);
//
//    ITensor t3(l1(2),l3(1),l2(1));
//
//    CHECK_EQUAL(t3.r(),3);
//    CHECK(hasindex(t3,l1));
//    CHECK(hasindex(t3,l2));
//    CHECK(hasindex(t3,l3));
//    CHECK_CLOSE(t3(l1(1),l3(1),l2(1)),0,1E-10);
//    CHECK_CLOSE(t3(l1(2),l3(1),l2(1)),1,1E-10);
//    CHECK_CLOSE(t3(l1(1),l3(2),l2(1)),0,1E-10);
//    CHECK_CLOSE(t3(l1(2),l3(2),l2(2)),0,1E-10);
//    CHECK_CLOSE(t3(l1(1),l3(1),l2(2)),0,1E-10);
//    CHECK_CLOSE(t3(l1(2),l3(1),l2(2)),0,1E-10);
//    CHECK_CLOSE(t3(l1(1),l3(2),l2(2)),0,1E-10);
//    CHECK_CLOSE(t3(l1(2),l3(2),l2(2)),0,1E-10);
//    CHECK_CLOSE(sumels(t3),1,1E-10);
//    CHECK_CLOSE(t3.norm(),1,1E-10);
//
//    ITensor t4(a1(1),l3(2),l2(1));
//
//    CHECK_EQUAL(t4.r(),3);
//    CHECK(hasindex(t4,a1));
//    CHECK(hasindex(t4,l2));
//    CHECK(hasindex(t4,l3));
//    CHECK_CLOSE(t4(l3(1),l2(1)),0,1E-10);
//    CHECK_CLOSE(t4(l3(1),l2(2)),0,1E-10);
//    CHECK_CLOSE(t4(l3(2),l2(1)),1,1E-10);
//    CHECK_CLOSE(t4(l3(2),l2(2)),0,1E-10);
//    CHECK_CLOSE(sumels(t4),1,1E-10);
//    CHECK_CLOSE(t4.norm(),1,1E-10);
//
//    ITensor r4(l1(1),l3(1),l2(2),l4(1));
//
//    CHECK_EQUAL(r4.r(),4);
//    CHECK(hasindex(r4,l1));
//    CHECK(hasindex(r4,l2));
//    CHECK(hasindex(r4,l3));
//    CHECK(hasindex(r4,l4));
//    CHECK_CLOSE(r4(l1(1),l3(1),l2(2),l4(1)),1,1E-10);
//    CHECK_CLOSE(sumels(r4),1,1E-10);
//    CHECK_CLOSE(r4.norm(),1,1E-10);
//
//    ITensor t8(l1(1),l2(2),l3(1),l4(2),l5(1),l6(2),l7(1),l8(2));
//
//    CHECK_EQUAL(t8.r(),8);
//    CHECK(hasindex(t8,l1));
//    CHECK(hasindex(t8,l2));
//    CHECK(hasindex(t8,l3));
//    CHECK(hasindex(t8,l4));
//    CHECK(hasindex(t8,l5));
//    CHECK(hasindex(t8,l6));
//    CHECK(hasindex(t8,l7));
//    CHECK(hasindex(t8,l8));
//
//    CHECK_CLOSE(t8(l1(1),l2(2),l3(1),l4(2),l5(1),l6(2),l7(1),l8(2)),1,1E-10);
//    CHECK_CLOSE(t8.norm(),1,1E-10);
//
//    }
//
//SECTION("MultiIndexConstructors")
//    {
//    const
//    IndexSet<Index> indices(a2,l3,l1,a4);
//
//    ITensor t1(indices);
//
//    CHECK_EQUAL(t1.r(),4);
//    CHECK(hasindex(t1,a2));
//    CHECK(hasindex(t1,l3));
//    CHECK(hasindex(t1,l1));
//    CHECK(hasindex(t1,a4));
//    CHECK_CLOSE(t1.norm(),0,1E-10);
//
//    Vector V(l1.m()*l3.m());
//    V.Randomize();
//
//    ITensor t2(indices,V);
//
//    CHECK_EQUAL(t2.r(),4);
//    CHECK(hasindex(t2,a2));
//    CHECK(hasindex(t2,l3));
//    CHECK(hasindex(t2,l1));
//    CHECK(hasindex(t2,a4));
//    CHECK_CLOSE(t2.norm(),Norm(V),1E-10);
//    CHECK_CLOSE(sumels(t2),V.sumels(),1E-10);
//    }
//
//SECTION("ITensorConstructors")
//{
//    Index clink("clink",4);
//    IndexSet<Index> indices1(l1,l2,clink);
//
//    Vector V(l1.m()*l2.m()*clink.m());
//    V.Randomize();
//
//    ITensor t1(indices1,V);
//
//    Real f = Global::random();
//
//    ITensor t2(t1);
//    t2 *= f;
//
//    IndexSet<Index> indices3(l1,l2,l3,l4);
//
//    ITensor t3(indices3,t2);
//
//    CHECK_EQUAL(4,t3.r());
//
//    for(int i = 1; i <= l1.m(); ++i)
//    for(int j = 1; j <= l2.m(); ++j)
//    {
//    CHECK_CLOSE(t1(l1(i),l2(j),clink(1))*f,t3(l1(i),l2(j),l3(1),l4(1)),1E-10);
//    CHECK_CLOSE(t1(l1(i),l2(j),clink(2))*f,t3(l1(i),l2(j),l3(2),l4(1)),1E-10);
//    CHECK_CLOSE(t1(l1(i),l2(j),clink(3))*f,t3(l1(i),l2(j),l3(1),l4(2)),1E-10);
//    CHECK_CLOSE(t1(l1(i),l2(j),clink(4))*f,t3(l1(i),l2(j),l3(2),l4(2)),1E-10);
//    }
//
//    Permutation P(NMAX+1);
//    P.setFromTo(2,4);
//    P.setFromTo(4,2);
//    //CHECK(P.check(4));
//
//    IndexSet<Index> indices5(l1,l4,l3,l2);
//
//    ITensor t4(t3);
//    Real f2 = Global::random();
//    t4 /= f2;
//    ITensor t5(indices5,t4,P);
//
//    CHECK_EQUAL(4,t5.r());
//
//    for(int i = 1; i <= l1.m(); ++i)
//    for(int j = 1; j <= l2.m(); ++j)
//    for(int k = 1; k <= l3.m(); ++k)
//    for(int l = 1; l <= l4.m(); ++l)
//    {
//    CHECK_CLOSE(t3(l1(i),l2(j),l3(k),l4(l))/f2,t5(l1(i),l2(j),l3(k),l4(l)),1E-10);
//    }
//
//}
//
//SECTION("Copy")
//{
//    IndexSet<Index> indices(a2,l3,l1,a4);
//
//    Vector V(l1.m()*l3.m());
//    V.Randomize();
//
//    ITensor t1(indices,V);
//
//    CHECK_EQUAL(t1.r(),4);
//    CHECK(hasindex(t1,a2));
//    CHECK(hasindex(t1,l3));
//    CHECK(hasindex(t1,l1));
//    CHECK(hasindex(t1,a4));
//    CHECK_CLOSE(t1.norm(),Norm(V),1E-10);
//    CHECK_CLOSE(sumels(t1),V.sumels(),1E-10);
//
//    //Use copy constructor
//    ITensor t2(t1);
//    t1 = ITensor(); //destroy t1
//
//    CHECK_EQUAL(t2.r(),4);
//    CHECK(hasindex(t2,a2));
//    CHECK(hasindex(t2,l3));
//    CHECK(hasindex(t2,l1));
//    CHECK(hasindex(t2,a4));
//    CHECK_CLOSE(t2.norm(),Norm(V),1E-10);
//    CHECK_CLOSE(sumels(t2),V.sumels(),1E-10);
//
//    //Use operator=
//    ITensor t3 = t2;
//    t2 = ITensor(); //destroy t2
//
//    CHECK_EQUAL(t3.r(),4);
//    CHECK(hasindex(t3,a2));
//    CHECK(hasindex(t3,l3));
//    CHECK(hasindex(t3,l1));
//    CHECK(hasindex(t3,a4));
//    CHECK_CLOSE(t3.norm(),Norm(V),1E-10);
//    CHECK_CLOSE(sumels(t3),V.sumels(),1E-10);
//}
//
//SECTION("ScalarMultiply")
//{
//    A *= -1;
//    CHECK_EQUAL(A(s1(1),s2(1)),-11);
//    CHECK_EQUAL(A(s1(1),s2(2)),-12);
//    CHECK_EQUAL(A(s1(2),s2(1)),-21);
//    CHECK_EQUAL(A(s1(2),s2(2)),-22);
//
//    Real f = Global::random();
//    A *= -f;
//    CHECK_CLOSE(A(s1(1),s2(1)),11*f,1E-10);
//    CHECK_CLOSE(A(s1(1),s2(2)),12*f,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(1)),21*f,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(2)),22*f,1E-10);
//
//    B /= f;
//    CHECK_CLOSE(B(s1(1),s2(1)),110/f,1E-10);
//    CHECK_CLOSE(B(s1(1),s2(2)),120/f,1E-10);
//    CHECK_CLOSE(B(s1(2),s2(1)),210/f,1E-10);
//    CHECK_CLOSE(B(s1(2),s2(2)),220/f,1E-10);
//}
//
///*
//SECTION("assignToVec")
//{
//    Vector V(l1.m()*l2.m()*l3.m());
//    V.Randomize();
//    Real f = -Global::random();
//
//    IndexSet<Index> indices(l1,l2,l3);
//
//    ITensor T(indices,V);
//
//    T *= f;
//
//    Vector U(T.indices().dim()); T.assignToVec(U);
//
//    CHECK_EQUAL(U.Length(),V.Length());
//
//    for(int j = 1; j < V.Length(); ++j)
//    { CHECK_CLOSE(U(j),V(j)*f,1E-10); }
//
//}
//*/
//
//SECTION("MapElems")
//    {
//    // class Functor and the function Func
//    // are defined in test.h
//
//    ITensor A1(A);
//    Functor f;
//    A1.mapElems(f);
//    for(int n1 = 1; n1 <= s1.m(); ++n1)
//    for(int n2 = 1; n2 <= s2.m(); ++n2)
//        {
//        CHECK_CLOSE( f( A(s1(n1),s2(n2)) ), A1(s1(n1),s2(n2)) ,1E-10);
//        }
//
//    ITensor A2(A);
//    Real (*pFunc)(Real) = Func;
//    A2.mapElems(*pFunc);
//    for(int n1 = 1; n1 <= s1.m(); ++n1)
//    for(int n2 = 1; n2 <= s2.m(); ++n2)
//        {
//        CHECK_CLOSE( Func( A(s1(n1),s2(n2)) ), A2(s1(n1),s2(n2)) ,1E-10);
//        }
//    }
//
///*
//SECTION("reshapeDat")
//    {
//    Permutation P;
//    P.fromTo(1,2);
//    P.fromTo(2,1);
//
//    Real f = -5;
//    A *= f;
//
//    A.reshapeDat(P);
//
//    CHECK_CLOSE(A(s1(1),s2(1)),11*f,1E-10);
//    CHECK_CLOSE(A(s1(1),s2(2)),21*f,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(1)),12*f,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(2)),22*f,1E-10);
//
//    }
//    */
//
///*
//SECTION("reshape")
//    {
//    //cout << "Begin: reshape -------------" << endl;
//    Permutation P;
//    P.fromTo(1,2);
//    P.fromTo(2,1);
//
//    Real f = -5;
//    A *= f;
//
//    CHECK_CLOSE(A(s1(1),s2(1)),11*f,1E-10);
//    CHECK_CLOSE(A(s1(1),s2(2)),12*f,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(1)),21*f,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(2)),22*f,1E-10);
//
//    ITensor cA(A);
//    cA.reshape(P);
//
//    //Check that indices are re-ordered...
//    CHECK(A.index(1) == cA.index(P.dest(1)));
//    CHECK(A.index(2) == cA.index(P.dest(2)));
//
//    //...but cA is equivalent to A apart from
//    //Index order:
//    CHECK_CLOSE(cA(s1(1),s2(1)),A(s1(1),s2(1)),1E-10);
//    CHECK_CLOSE(cA(s1(1),s2(2)),A(s1(1),s2(2)),1E-10);
//    CHECK_CLOSE(cA(s1(2),s2(1)),A(s1(2),s2(1)),1E-10);
//    CHECK_CLOSE(cA(s1(2),s2(2)),A(s1(2),s2(2)),1E-10);
//
//    //cout << "End: reshape ---------------" << endl;
//    }
//    */
//
//SECTION("SumDifference")
//{
//    Vector V(mixed_inds_dim),W(mixed_inds_dim);
//    V.Randomize();
//    W.Randomize();
//
//    ITensor v(mixed_inds,V), w(mixed_inds,W);
//
//    Real f1 = -Global::random(), f2 = 0.1*f1;
//
//    ITensor r = f1*v + w/f2; 
//    VectorRef R = r.assignToVec();
//    for(int j = 1; j < R.Length(); ++j)
//    { CHECK_CLOSE(R(j),f1*V(j)+W(j)/f2,1E-10); }
//
//    ITensor d(v); d -= w;
//    VectorRef D = d.assignToVec();
//    for(int j = 1; j < D.Length(); ++j)
//    { CHECK_CLOSE(D(j),V(j)-W(j),1E-10); }
//
//    Vector YY(mixed_inds_dim),ZZ(mixed_inds_dim);
//    YY.Randomize();
//    ZZ.Randomize();
//
//    f1 = 1; f2 = 1;
//    ITensor yy(mixed_inds,YY), zz(mixed_inds,ZZ);
//    r = f1*yy + f2*zz;
//    for(int j1 = 1; j1 <= 2; ++j1)
//    for(int j2 = 1; j2 <= 2; ++j2)
//    for(int k3 = 1; k3 <= 3; ++k3)
//    for(int j4 = 1; j4 <= 2; ++j4)
//    { 
//        CHECK_CLOSE(r(l1(j1),l2(j2),b3(k3),l4(j4)),
//                    f1*yy(l1(j1),l2(j2),b3(k3),l4(j4))
//                   +f2*zz(l1(j1),l2(j2),b3(k3),l4(j4))
//                   ,1E-10); 
//    }
//
//}
//
//SECTION("ContractingProduct")
//    {
//
//    //Check for rank 0 ITensors
//    {
//    Real f = Global::random();
//    ITensor rZ(f), T(b2,a1,b4);
//    T.randomize();
//
//    ITensor res = rZ * T;
//
//    CHECK_EQUAL(rZ.r(),0);
//    CHECK_EQUAL(res.r(),3);
//
//    for(int j2 = 1; j2 <= 2; ++j2)
//    for(int j4 = 1; j4 <= 4; ++j4)
//        {
//        Real val = f * T(b2(j2),a1(1),b4(j4));
//        CHECK_CLOSE(res(b2(j2),a1(1),b4(j4)),val,1E-10);
//        }
//    }
//
//    //More general case
//    ITensor L(b4,a1,b3,a2,b2), R(b5,a1,b4,b2,b3);
//
//    L.randomize(); R.randomize();
//
//    Real fL = Global::random(), fR = Global::random();
//    ITensor Lf = L * fL;
//    ITensor Rf = R * fR;
//
//    ITensor res1 = Lf*Rf;
//
//    CHECK(hasindex(res1,b5));
//    CHECK(hasindex(res1,a2));
//    CHECK(!hasindex(res1,a1));
//    CHECK(!hasindex(res1,b2));
//    CHECK(!hasindex(res1,b3));
//    CHECK(!hasindex(res1,b4));
//
//    CHECK_EQUAL(res1.r(),2);
//
//    for(int j5 = 1; j5 <= b5.m(); ++j5)
//        {
//        Real val = 0;
//        for(int j2 = 1; j2 <= 2; ++j2)
//        for(int j3 = 1; j3 <= 3; ++j3)
//        for(int j4 = 1; j4 <= 4; ++j4)
//            {
//            val += L(b2(j2),a1(1),b3(j3),b4(j4))*fL * R(b5(j5),a1(1),b3(j3),b2(j2),b4(j4))*fR;
//            }
//        CHECK_CLOSE(res1(a2(1),b5(j5)),val,1E-10);
//        }
//
//    ITensor res2 = R*L;
//
//    CHECK(hasindex(res2,b5));
//    CHECK(hasindex(res2,a2));
//    CHECK(!hasindex(res2,a1));
//    CHECK(!hasindex(res2,b2));
//    CHECK(!hasindex(res2,b3));
//    CHECK(!hasindex(res2,b4));
//
//    CHECK_EQUAL(res2.r(),2);
//
//    for(int j5 = 1; j5 <= b5.m(); ++j5)
//        {
//        Real val = 0;
//        for(int j2 = 1; j2 <= 2; ++j2)
//        for(int j3 = 1; j3 <= 3; ++j3)
//        for(int j4 = 1; j4 <= 4; ++j4)
//            {
//            val += L(b2(j2),a1(1),b3(j3),b4(j4)) * R(b5(j5),a1(1),b3(j3),b2(j2),b4(j4));
//            }
//        CHECK_CLOSE(res2(a2(1),b5(j5)),val,1E-10);
//        }
//
//    ITensor Q(a1,b4,a2,b2), P(a2,a3,a1);
//
//    Q.randomize(); P.randomize();
//
//    Real fQ = Global::random(), fP = Global::random();
//    ITensor Qf = Q * fQ;
//    ITensor Pf = P * fP;
//
//    ITensor res3 = Qf*Pf;
//
//    CHECK(hasindex(res3,b4));
//    CHECK(hasindex(res3,b2));
//    CHECK(hasindex(res3,a3));
//    CHECK(!hasindex(res3,a1));
//    CHECK(!hasindex(res3,a2));
//
//    CHECK_EQUAL(res3.r(),3);
//
//    for(int j2 = 1; j2 <= b2.m(); ++j2)
//    for(int j4 = 1; j4 <= b4.m(); ++j4)
//        {
//        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
//        CHECK_CLOSE(res3(b4(j4),b2(j2)),val,1E-10);
//        }
//
//    ITensor res4 = Pf*Qf;
//
//    CHECK(hasindex(res4,b4));
//    CHECK(hasindex(res4,b2));
//    CHECK(hasindex(res4,a3));
//    CHECK(!hasindex(res4,a1));
//    CHECK(!hasindex(res4,a2));
//
//    CHECK_EQUAL(res4.r(),3);
//
//    for(int j2 = 1; j2 <= 2; ++j2)
//    for(int j4 = 1; j4 <= 4; ++j4)
//        {
//        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
//        CHECK_CLOSE(res4(b4(j4),b2(j2)),val,1E-10);
//        }
//
//
//    ITensor psi(a1,a2,a3), mpoh(l2,a1,prime(a1),a2,prime(a2));
//    psi.randomize(); mpoh.randomize();
//
//    ITensor Hpsi = mpoh * psi;
//
//    CHECK_EQUAL(Hpsi.r(),4);
//    CHECK(hasindex(Hpsi,l2));
//    CHECK(hasindex(Hpsi,prime(a1)));
//    CHECK(hasindex(Hpsi,prime(a2)));
//    CHECK(hasindex(Hpsi,a3));
//    CHECK(!hasindex(Hpsi,a1));
//    CHECK(!hasindex(Hpsi,a2));
//    }
//
//SECTION("NonContractingProduct")
//    {
//    ITensor L(b2,a1,b3,b4), R(a1,b3,a2,b5,b4);
//
//    L.randomize(); R.randomize();
//
//    Real fL = Global::random(), fR = Global::random();
//    ITensor Lf = L * fL;
//    ITensor Rf = R * fR;
//
//    ITensor res1 = Lf / Rf;
//
//    CHECK(hasindex(res1,b2));
//    CHECK(hasindex(res1,a2));
//    CHECK(hasindex(res1,b5));
//    CHECK(hasindex(res1,a1));
//    CHECK(hasindex(res1,b3));
//    CHECK(hasindex(res1,b4));
//
//    CHECK_EQUAL(res1.r(),6);
//
//    for(int j2 = 1; j2 <= 2; ++j2)
//    for(int j3 = 1; j3 <= 3; ++j3)
//    for(int j4 = 1; j4 <= 4; ++j4)
//    for(int j5 = 1; j5 <= 5; ++j5)
//    {
//        Real val = L(b2(j2),b3(j3),b4(j4))*fL * R(b3(j3),b5(j5),b4(j4))*fR;
//        CHECK_CLOSE(res1(b2(j2),b3(j3),b4(j4),b5(j5)),val,1E-10);
//    }
//
//    ITensor res2 = R/L;
//
//    CHECK(hasindex(res2,b2));
//    CHECK(hasindex(res2,a2));
//    CHECK(hasindex(res2,b5));
//    CHECK(hasindex(res2,a1));
//    CHECK(hasindex(res2,b3));
//    CHECK(hasindex(res2,b4));
//
//    CHECK_EQUAL(res2.r(),6);
//
//    for(int j2 = 1; j2 <= 2; ++j2)
//    for(int j3 = 1; j3 <= 3; ++j3)
//    for(int j4 = 1; j4 <= 4; ++j4)
//    for(int j5 = 1; j5 <= 5; ++j5)
//    {
//        Real val = L(b2(j2),a1(1),b3(j3),b4(j4)) * R(a1(1),b3(j3),a2(1),b5(j5),b4(j4));
//        CHECK_CLOSE(res2(b2(j2),b3(j3),b4(j4),b5(j5)),val,1E-10);
//    }
//
//    ITensor Q(a1,b4,a2,b2), P(a2,a3,a1);
//
//    Q.randomize(); P.randomize();
//
//    Real fQ = Global::random(), fP = Global::random();
//    ITensor Qf = Q * fQ;
//    ITensor Pf = P * fP;
//
//    ITensor res3 = Qf/Pf;
//
//    CHECK(hasindex(res3,b4));
//    CHECK(hasindex(res3,b2));
//    CHECK(hasindex(res3,a3));
//    CHECK(hasindex(res3,a1));
//    CHECK(hasindex(res3,a2));
//
//    CHECK_EQUAL(res3.r(),5);
//
//    for(int j2 = 1; j2 <= 2; ++j2)
//    for(int j4 = 1; j4 <= 4; ++j4)
//    {
//        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
//        CHECK_CLOSE(res3(b4(j4),b2(j2)),val,1E-10);
//    }
//
//    ITensor res4 = Pf/Qf;
//
//    CHECK(hasindex(res4,b4));
//    CHECK(hasindex(res4,b2));
//    CHECK(hasindex(res4,a3));
//    CHECK(hasindex(res4,a1));
//    CHECK(hasindex(res4,a2));
//
//    CHECK_EQUAL(res4.r(),5);
//
//    for(int j2 = 1; j2 <= 2; ++j2)
//    for(int j4 = 1; j4 <= 4; ++j4)
//    {
//        Real val = Q(a1(1),b4(j4),a2(1),b2(j2))*fQ * P(a2(1),a3(1),a1(1))*fP;
//        CHECK_CLOSE(res4(b4(j4),b2(j2)),val,1E-10);
//    }
//
//
//    ITensor psi(a1,a2,a3), mpoh(l2,a1,prime(a1),a2,prime(a2));
//    psi.randomize(); mpoh.randomize();
//
//    ITensor Hpsi = mpoh / psi;
//
//    CHECK_EQUAL(Hpsi.r(),6);
//    CHECK(hasindex(Hpsi,l2));
//    CHECK(hasindex(Hpsi,a1));
//    CHECK(hasindex(Hpsi,a2));
//    CHECK(hasindex(Hpsi,prime(a1)));
//    CHECK(hasindex(Hpsi,prime(a2)));
//    CHECK(hasindex(Hpsi,a3));
//
//    for(int j2 = 1; j2 <= 2; ++j2)
//        { CHECK_CLOSE(Hpsi(l2(j2)),psi(a1(1),a2(1),a3(1))*mpoh(l2(j2)),1E-10); }
//    }
//
//SECTION("ComplexNonContractingProduct")
//    {
//    ITensor Lr(b2,b3,b4), Li(b2,b3,b4),
//            Rr(b3,a2,b5,b4), Ri(b3,a2,b5,b4);
//
//    Lr.randomize(); 
//    Li.randomize(); 
//    Rr.randomize();
//    Ri.randomize();
//
//    ITensor L = Complex_1*Lr + Complex_i*Li;
//    ITensor R = Complex_1*Rr + Complex_i*Ri;
//
//    ITensor res1 = L / R;
//
//    CHECK(hasindex(res1,b2));
//    CHECK(hasindex(res1,a2));
//    CHECK(hasindex(res1,b5));
//    CHECK(hasindex(res1,b3));
//    CHECK(hasindex(res1,b4));
//
//    CHECK_EQUAL(res1.r(),5);
//
//    ITensor resR(realPart(res1)),
//            resI(imagPart(res1));
//
//    ITensor rdiff = resR-(Lr/Rr-Li/Ri);
//    ITensor idiff = resI-(Lr/Ri+Li/Rr);
//
//    CHECK(rdiff.norm() < 1E-12);
//    CHECK(idiff.norm() < 1E-12);
//
//    }
//
//SECTION("TieIndices")
//    {
//
//    Index t("tied",2);
//
//    ITensor dX(X);
//    dX.tieIndices(s1,s2,t);
//
//    CHECK_CLOSE(dX.norm(),0,1E-5);
//    CHECK_EQUAL(dX.r(),1);
//    CHECK(hasindex(dX,t));
//
//    ITensor dZ(Z);
//    dZ.tieIndices(s1,s2,t);
//    CHECK_CLOSE(dZ(t(1)),+1,1E-5);
//    CHECK_CLOSE(dZ(t(2)),-1,1E-5);
//
//    {
//    ITensor T(l1,l2,a1,s2,s1);
//    T.randomize();
//
//    ITensor TT(T);
//    TT.tieIndices(l2,l1,s1,l2);
//
//    CHECK_EQUAL(TT.r(),3);
//
//    for(int j = 1; j <= 2; ++j)
//    for(int k = 1; k <= 2; ++k)
//        {
//        CHECK_CLOSE(T(l1(j),l2(j),a1(1),s2(k),s1(j)),TT(l2(j),s2(k),a1(1)),1E-5);
//        }
//    }
//
//    //Try tying m==1 inds
//    {
//    ITensor T(l1,a2,a1,s2,a3);
//    T.randomize();
//
//    ITensor TT(T);
//    TT.tieIndices(a1,a3,a2,a1);
//
//    CHECK_EQUAL(TT.r(),3);
//
//    for(int j = 1; j <= 2; ++j)
//    for(int k = 1; k <= 2; ++k)
//        {
//        CHECK_CLOSE(T(l1(j),s2(k)),TT(l1(j),s2(k)),1E-5);
//        }
//    }
//
//
//
//    } //TieIndices
//
//SECTION("ComplexTieIndices")
//    {
//    ITensor Tr(l1,l2,a1,s2,s1),
//            Ti(l1,l2,a1,s2,s1);
//    Tr.randomize();
//    Ti.randomize();
//
//    ITensor T = Tr + Complex_i*Ti;
//
//    ITensor TT(T);
//    TT.tieIndices(l2,l1,s1,l2);
//
//    CHECK_EQUAL(TT.r(),3);
//
//    ITensor TTr(realPart(TT)),
//            TTi(imagPart(TT));
//
//    for(int j = 1; j <= 2; ++j)
//    for(int k = 1; k <= 2; ++k)
//        {
//        CHECK_CLOSE(Tr(l1(j),l2(j),a1(1),s2(k),s1(j)),TTr(l2(j),s2(k),a1(1)),1E-5);
//        CHECK_CLOSE(Ti(l1(j),l2(j),a1(1),s2(k),s1(j)),TTi(l2(j),s2(k),a1(1)),1E-5);
//        }
//    }
//
//SECTION("Trace")
//    {
//
//    ITensor A(b2,a1,b3,b5,prime(b3));
//    A.randomize();
//    Real f = -Global::random();
//    A *= f;
//
//    ITensor At = trace(A,b3,prime(b3));
//
//    for(int j2 = 1; j2 <= b2.m(); ++j2)
//    for(int j5 = 1; j5 <= b5.m(); ++j5)
//        {
//        Real val = 0;
//        for(int j3 = 1; j3 <= b3.m(); ++j3)
//            {
//            val += A(b2(j2),a1(1),b3(j3),b5(j5),prime(b3)(j3));
//            }
//        CHECK_CLOSE(val,At(b2(j2),a1(1),b5(j5)),1E-10);
//        }
//
//    ITensor MM(b5,prime(b5));
//    MM.randomize();
//    MM *= -2.34;
//
//    Real tr = trace(MM);
//
//    Real check_tr = 0;
//    for(int j5 = 1; j5 <= b5.m(); ++j5)
//        {
//        check_tr += MM(b5(j5),prime(b5)(j5));
//        }
//    CHECK_CLOSE(tr,check_tr,1E-10);
//
//    }
//
//SECTION("fromMatrix11")
//    {
//    Matrix M22(s1.m(),s2.m());
//
//    M22(1,1) = -0.3; M22(1,2) = 110;
//    M22(2,1) = -1.7; M22(1,2) = 5;
//
//    ITensor T(s1,s2);
//    //T should be overwritten so check
//    //that scalar mult has no effect
//    T *= -5; 
//
//    T.fromMatrix11(s1,s2,M22);
//
//    CHECK_CLOSE(T(s1(1),s2(1)),M22(1,1),1E-10);
//    CHECK_CLOSE(T(s1(1),s2(2)),M22(1,2),1E-10);
//    CHECK_CLOSE(T(s1(2),s2(1)),M22(2,1),1E-10);
//    CHECK_CLOSE(T(s1(2),s2(2)),M22(2,2),1E-10);
//
//    ITensor U(T);
//
//    U.fromMatrix11(s2,s1,M22);
//
//    CHECK_CLOSE(T(s1(1),s2(1)),M22(1,1),1E-10);
//    CHECK_CLOSE(T(s1(1),s2(2)),M22(1,2),1E-10);
//    CHECK_CLOSE(T(s1(2),s2(1)),M22(2,1),1E-10);
//    CHECK_CLOSE(T(s1(2),s2(2)),M22(2,2),1E-10);
//
//    CHECK_CLOSE(U(s2(1),s1(1)),M22(1,1),1E-10);
//    CHECK_CLOSE(U(s2(1),s1(2)),M22(1,2),1E-10);
//    CHECK_CLOSE(U(s2(2),s1(1)),M22(2,1),1E-10);
//    CHECK_CLOSE(U(s2(2),s1(2)),M22(2,2),1E-10);
//
//    Matrix M12(a1.m(),s2.m());
//    M12(1,1) = 37; M12(1,2) = -2;
//
//    ITensor P(a1,s2);
//    P *= -4;
//
//    P.fromMatrix11(a1,s2,M12);
//
//    CHECK_CLOSE(P(a1(1),s2(1)),M12(1,1),1E-10);
//    CHECK_CLOSE(P(a1(1),s2(2)),M12(1,2),1E-10);
//
//    P.fromMatrix11(s2,a1,M12.t());
//
//    CHECK_CLOSE(P(s2(1),a1(1)),M12(1,1),1E-10);
//    CHECK_CLOSE(P(s2(2),a1(1)),M12(1,2),1E-10);
//    }
//
//SECTION("ToFromMatrix11")
//    {
//    Matrix M(s1.m(),s2.m());    
//
//    Real f = -Global::random();
//
//    A *= f;
//
//    A.toMatrix11(s2,s1,M);
//
//    CHECK_CLOSE(M(1,1),11*f,1E-10);
//    CHECK_CLOSE(M(2,1),12*f,1E-10);
//    CHECK_CLOSE(M(1,2),21*f,1E-10);
//    CHECK_CLOSE(M(2,2),22*f,1E-10);
//
//    A.toMatrix11(s1,s2,M);
//
//    CHECK_CLOSE(M(1,1),11*f,1E-10);
//    CHECK_CLOSE(M(1,2),12*f,1E-10);
//    CHECK_CLOSE(M(2,1),21*f,1E-10);
//    CHECK_CLOSE(M(2,2),22*f,1E-10);
//
//    A.toMatrix11NoScale(s2,s1,M);
//
//    CHECK_CLOSE(M(1,1),11,1E-10);
//    CHECK_CLOSE(M(2,1),12,1E-10);
//    CHECK_CLOSE(M(1,2),21,1E-10);
//    CHECK_CLOSE(M(2,2),22,1E-10);
//
//    A *= -40;
//    A.fromMatrix11(s2,s1,M);
//
//    CHECK_CLOSE(A(s1(1),s2(1)),11,1E-10);
//    CHECK_CLOSE(A(s1(1),s2(2)),12,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(1)),21,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(2)),22,1E-10);
//
//    A.fromMatrix11(s1,s2,M);
//
//    CHECK_CLOSE(A(s1(1),s2(1)),11,1E-10);
//    CHECK_CLOSE(A(s1(1),s2(2)),21,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(1)),12,1E-10);
//    CHECK_CLOSE(A(s1(2),s2(2)),22,1E-10);
//
//
//    Vector V(4);
//    V(1) = 3.14; V(2) = 2.718; V(3) = -1; V(4) = 0;
//    Index link("link",4);
//
//    ITensor T(link,a1);
//    T(link(1),a1(1)) = V(1);
//    T(link(2),a1(1)) = V(2);
//    T(link(3),a1(1)) = V(3);
//    T(link(4),a1(1)) = V(4);
//
//    Matrix M41(4,1), M14(1,4);
//    
//    T.toMatrix11(link,a1,M41);
//
//    CHECK_CLOSE(M41(1,1),V(1),1E-10);
//    CHECK_CLOSE(M41(2,1),V(2),1E-10);
//    CHECK_CLOSE(M41(3,1),V(3),1E-10);
//    CHECK_CLOSE(M41(4,1),V(4),1E-10);
//     
//    T.toMatrix11(a1,link,M14);
//
//    CHECK_CLOSE(M14(1,1),V(1),1E-10);
//    CHECK_CLOSE(M14(1,2),V(2),1E-10);
//    CHECK_CLOSE(M14(1,3),V(3),1E-10);
//    CHECK_CLOSE(M14(1,4),V(4),1E-10);
//
//    }
//
//SECTION("ToFromMatrix22")
//    {
//    Index i1("i1",3),
//          i2("i2",4),
//          i3("i3",2),
//          i4("i4",4);
//
//    ITensor T(i1,i2,i3,i4);
//    T.randomize();
//    T *= -1.23235;
//
//    Matrix M;
//    T.toMatrix22(i2,i1,i4,i3,M);
//    ITensor V;
//    V.fromMatrix22(i2,i1,i4,i3,M);
//
//    CHECK((T-V).norm() < 1E-12);
//    }
//
///*
//SECTION("SymmetricDiag11")
//    {
//    ITensor T(s1,prime(s1));
//    commaInit(T,s1,prime(s1)) << 1, 2,
//                                  2, 1;
//
//    T *= -2;
//    Index mid;
//    ITensor D,U;
//    T.symmetricDiag11(s1,D,U,mid);
//    ITensor UD(U);
//    UD.prime(s1);
//    UD /= D;
//    ITensor rT = UD*U;
//    ITensor diff(UD*U - T);
//    CHECK(diff.norm() < 1E-10);
//
//    //Construct a random, symmetric ITensor
//    const int qs = 50;
//    Index q("q",qs);
//    Matrix QQ(qs,qs);
//    QQ.Randomize();
//    QQ += QQ.t();
//    ITensor Q(q,prime(q),QQ);
//    Q *= -2;
//
//    //Diagonalize and check the factorization
//    Q.symmetricDiag11(q,D,U,mid);
//    UD =U;
//    UD.prime(q);
//    UD /= D;
//    diff = UD*U - Q;
//    CHECK(diff.norm() < 1E-10);
//    }
//    */
//
//SECTION("CommaAssignment")
//    {
//    ITensor VV(s1);
//    VV.randomize();
//    VV *= -1;
//    commaInit(VV,s1) << 1, 2;
//    CHECK_EQUAL(VV(s1(1)),1);
//    CHECK_EQUAL(VV(s1(2)),2);
//
//    ITensor ZZ(s1,s2);
//    commaInit(ZZ,s1,s2) << 1, 0, 
//                           0, -1;
//    CHECK_EQUAL(ZZ(s1(1),s2(1)),1);
//    CHECK_EQUAL(ZZ(s1(2),s2(1)),0);
//    CHECK_EQUAL(ZZ(s1(1),s2(2)),0);
//    CHECK_EQUAL(ZZ(s1(2),s2(2)),-1);
//
//    ITensor XX(s1,s2);
//    XX(s1(2),s2(1)) = 5;
//    XX *= 3;
//    commaInit(XX,s1,s2) << 0, 1, 
//                           1, 0;
//    CHECK_EQUAL(XX(s1(1),s2(1)),0);
//    CHECK_EQUAL(XX(s1(2),s2(1)),1);
//    CHECK_EQUAL(XX(s1(1),s2(2)),1);
//    CHECK_EQUAL(XX(s1(2),s2(2)),0);
//
//    ITensor AA(s1,s2);
//    AA.randomize();
//    AA *= -Global::random();
//    commaInit(AA,s1,s2) << 11, 12, 
//                           21, 22;
//    CHECK_EQUAL(AA(s1(1),s2(1)),11);
//    CHECK_EQUAL(AA(s1(1),s2(2)),12);
//    CHECK_EQUAL(AA(s1(2),s2(1)),21);
//    CHECK_EQUAL(AA(s1(2),s2(2)),22);
//
//    ITensor T(s1,s2,s3);
//    T.randomize();
//    T *= -Global::random();
//    commaInit(T,s1,s2,s3) << 111, 112, 
//                             121, 122,
//                             211, 212,
//                             221, 222;
//    CHECK_EQUAL(T(s1(1),s2(1),s3(1)),111);
//    CHECK_EQUAL(T(s1(1),s2(1),s3(2)),112);
//    CHECK_EQUAL(T(s1(1),s2(2),s3(1)),121);
//    CHECK_EQUAL(T(s1(1),s2(2),s3(2)),122);
//    CHECK_EQUAL(T(s1(2),s2(1),s3(1)),211);
//    CHECK_EQUAL(T(s1(2),s2(1),s3(2)),212);
//    CHECK_EQUAL(T(s1(2),s2(2),s3(1)),221);
//    CHECK_EQUAL(T(s1(2),s2(2),s3(2)),222);
//    }
//
//SECTION("RealImagPart")
//    {
//    const Real f1 = 2.124,
//               f2 = 1.113;
//    ITensor ZiX = f1*Complex_1*Z + f2*Complex_i*X;
//    ITensor R(realPart(ZiX)),
//            I(imagPart(ZiX));
//    R -= f1*Z;
//    I -= f2*X;
//    CHECK_CLOSE(R.norm(),0,1E-5);
//    CHECK_CLOSE(I.norm(),0,1E-5);
//
//    //Test hc:
//    
//    ZiX.dag();
//    R = realPart(ZiX);
//    I = imagPart(ZiX);
//    R -= f1*Z;
//    I += f2*X;
//    CHECK_CLOSE(R.norm(),0,1E-5);
//    CHECK_CLOSE(I.norm(),0,1E-5);
//    }
//
//SECTION("SwapPrimeTest")
//    {
//    ITensor T(s1,prime(s1));
//    commaInit(T,s1,prime(s1)) << 11, 12,
//                                  21, 22;
//
//    CHECK_EQUAL(T(s1(1),prime(s1)(1)),11);
//    CHECK_EQUAL(T(s1(2),prime(s1)(1)),21);
//    CHECK_EQUAL(T(s1(1),prime(s1)(2)),12);
//    CHECK_EQUAL(T(s1(2),prime(s1)(2)),22);
//
//    T = swapPrime(T,0,1);
//
//    CHECK_EQUAL(T(prime(s1)(1),s1(1)),11);
//    CHECK_EQUAL(T(prime(s1)(2),s1(1)),21);
//    CHECK_EQUAL(T(prime(s1)(1),s1(2)),12);
//    CHECK_EQUAL(T(prime(s1)(2),s1(2)),22);
//    }
//
//SECTION("NoprimeTest")
//    {
//    ITensor T(s1,prime(s1));
//
//    //Check that T.noprime()
//    //throws an exception since it would
//    //lead to duplicate indices
//    CHECK_THROWS_AS(T.noprime(),ITError);
//    }
//
//SECTION("NormTest")
//    {
//    A.randomize();
//    CHECK_CLOSE(A.norm(),sqrt((A*A).toReal()),1E-5);
//
//    ITensor C = Complex_1*A+Complex_i*B;
//
//    CHECK_CLOSE(C.norm(),sqrt(realPart(dag(C)*C).toReal()),1E-5);
//    }
//
//SECTION("CR_ComplexAddition")
//    {
//    const Real f1 = 1.234,
//               f2 = 2.456;
//    ITensor iZX = f1*Complex_i*Z + f2*Complex_1*X;
//    ITensor R(realPart(iZX)),
//            I(imagPart(iZX));
//    R -= f2*X;
//    I -= f1*Z;
//    CHECK_CLOSE(R.norm(),0,1E-5);
//    CHECK_CLOSE(I.norm(),0,1E-5);
//    }
//
//SECTION("CC_ComplexAddition")
//    {
//    const Real f1 = 1.234,
//               f2 = 2.456;
//    ITensor iZiX = f1*Complex_i*Z + f2*Complex_i*X;
//    ITensor R(realPart(iZiX)),
//            I(imagPart(iZiX));
//    I -= f1*Z+f2*X;
//    CHECK_CLOSE(R.norm(),0,1E-5);
//    CHECK_CLOSE(I.norm(),0,1E-5);
//    }
//
//SECTION("ComplexScalar")
//    {
//    ITensor A(b4,s1),
//            B(b4,s1);
//    A.randomize();
//    B.randomize();
//    
//    const Real f1 = 2.1324,
//               f2 = -5.2235;
//
//    ITensor T1 = Complex(f1,0)*A;
//
//    CHECK((realPart(T1)-(f1*A)).norm() < 1E-12);
//    CHECK(imagPart(T1).norm() < 1E-12);
//
//    ITensor T2 = Complex(0,f2)*A;
//
//    CHECK(realPart(T2).norm() < 1E-12);
//    CHECK((imagPart(T2)-f2*A).norm() < 1E-12);
//
//    ITensor T3 = Complex(f1,f2)*A;
//    CHECK((realPart(T3)-f1*A).norm() < 1E-12);
//    CHECK((imagPart(T3)-f2*A).norm() < 1E-12);
//
//    ITensor T4 = Complex(f2,f1)*A;
//    CHECK((realPart(T4)-f2*A).norm() < 1E-12);
//    CHECK((imagPart(T4)-f1*A).norm() < 1E-12);
//
//    ITensor T5 = A+Complex_i*B;
//    T5 *= Complex(f1,f2);
//    CHECK((realPart(T5)-(f1*A-f2*B)).norm() < 1E-12);
//    CHECK((imagPart(T5)-(f2*A+f1*B)).norm() < 1E-12);
//
//    ITensor T6 = A+Complex_i*B;
//    T6 *= Complex(f2,f1);
//    CHECK((realPart(T6)-(f2*A-f1*B)).norm() < 1E-12);
//    CHECK((imagPart(T6)-(f1*A+f2*B)).norm() < 1E-12);
//    }
//
//SECTION("CommonIndex")
//    {
//    ITensor T1(s1,s2,l1,l2),
//            T2(s1,l3),
//            T3(s3,l4);
//
//    Index c = commonIndex(T1,T3);
//    CHECK(!c);
//
//    c = commonIndex(T2,T3);
//    CHECK(!c);
//
//    CHECK(commonIndex(T1,T2) == s1);
//    CHECK(commonIndex(T1,T2,Site) == s1);
//    }
//
//SECTION("DiagITensorBasicContraction")
//    {
//    Vector v(3);
//    v(1) = -0.8;
//    v(2) = 1.7;
//    v(3) = 4.9;
//
//    Vector vb(2);
//    vb(1) = 1;
//    vb(2) = -1;
//
//    const Real f1 = Global::random(),
//               f2 = Global::random();
//
//    ITensor op1(s1,prime(s1),f1),
//            op2(s1,prime(s1),f2),
//            opa(s1,a1,3.1),
//            psi(s1,l1,-1),
//            opb(s1,b2,vb),
//            r1(s1,prime(s1,2)),
//            r2(s1,prime(s1,2));
//
//    r1.randomize();
//    r2.randomize();
//
//    CHECK(op1.type() == ITensor::Diag);
//    CHECK(op2.type() == ITensor::Diag);
//    CHECK(opa.type() == ITensor::Diag);
//    CHECK(opb.type() == ITensor::Diag);
//    CHECK(psi.type() == ITensor::Diag);
//
//    CHECK(r1.type() == ITensor::Dense);
//    CHECK(r2.type() == ITensor::Dense);
//    
//    ITensor op3 = op1 + Complex_i*op2;
//
//    CHECK(op3.type() == ITensor::Diag);
//
//    ITensor res1 = op1*r1;
//    CHECK(res1.type() == ITensor::Dense);
//    res1.mapprime(1,0);
//    ITensor diff1 = res1-f1*r1;
//    CHECK(diff1.norm() < 1E-10);
//
//    ITensor res2 = r1*op1;
//    CHECK(res2.type() == ITensor::Dense);
//    res2.mapprime(1,0);
//    ITensor diff2 = res2-f1*r1;
//    CHECK(diff2.norm() < 1E-10);
//
//    ITensor rc = r1+Complex_i*r2;
//
//    ITensor res3 = rc*op1;
//    res3.mapprime(1,0);
//    CHECK(res3.type() == ITensor::Dense);
//    CHECK(res3.isComplex());
//    ITensor diff3 = res3-f1*rc;
//    CHECK(diff3.norm() < 1E-10);
//
//    ITensor res4 = op1*rc;
//    res4.mapprime(1,0);
//    CHECK(res4.type() == ITensor::Dense);
//    CHECK(res4.isComplex());
//    ITensor diff4 = res4-f1*rc;
//    CHECK(diff4.norm() < 1E-10);
//
//    ITensor res5 = rc*op3;
//    CHECK(res5.type() == ITensor::Dense);
//    CHECK(res5.isComplex());
//    ITensor rres5(realPart(res5)),
//            ires5(imagPart(res5));
//    ITensor rdiff5 = rres5-(r1*op1-r2*op2),
//            idiff5 = ires5-(r1*op2+r2*op1);
//    CHECK(rdiff5.norm() < 1E-10);
//    CHECK(idiff5.norm() < 1E-10);
//    }
//
//SECTION("DiagMethod")
//    {
//    ITensor t1(b3,b4);
//    t1.randomize();
//    t1 *= -8.232244;
//    Vector d1 = t1.diag();
//    for(int i = 1; i <= minM(t1.indices()); ++i)
//        {
//        CHECK(fabs(d1(i)-t1(b3(i),b4(i))) < 1E-12);
//        }
//
//    Vector v(4);
//    v(1) = -2.2442;
//    v(2) = 1.34834;
//    v(3) = 0.0;
//    v(4) = 8.38457;
//    ITensor t2(prime(b4),b4,v);
//    CHECK(t2.type() == ITensor::Diag);
//    CHECK(Norm(v-t2.diag()) < 1E-12);
//    }

}
