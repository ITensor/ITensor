#include "test.h"
#include "iqtensor.h"

using namespace itensor;
using namespace std;

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


TEST_CASE("IQTensorTest")
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

IQIndex S1,S2,L1,L2;

IQTensor phi,A,B,C,D;

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

phi = IQTensor(S1(1),S2(1),L2(3));
phi.randomize();

A = IQTensor(L1(1),S1(1),L2(4),S2(2));
A.randomize();

B = IQTensor(L1(1),L2(3));
B.randomize();

C = IQTensor(dag(L1)(5),prime(L1)(5));
C.randomize();

D = IQTensor(dag(L1)(3),S1(1),prime(L1)(3),prime(L1,2)(5));
D.randomize();

SECTION("Boolean")
    {
    IQTensor t1;

    CHECK(!t1);

    CHECK(A);
    CHECK(B);
    CHECK(C);
    CHECK(D);
    }

SECTION("Constructors")
    {
    Real f = Global::random();
    IQTensor rZ(f);

    CHECK_EQUAL(rZ.r(),0);
    CHECK_CLOSE(rZ.norm(),f,1E-10);
    }

SECTION("Contracting Product")
    {
    IQTensor res = A * dag(B);

    CHECK(hasindex(res,S1));
    CHECK(hasindex(res,S2));
    CHECK(!hasindex(res,L1));
    CHECK(!hasindex(res,L2));

    for(int k1 = 1; k1 <= S1.m(); ++k1)
    for(int k2 = 1; k2 <= S2.m(); ++k2)
        {
        Real val = 0;
        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
            {
            val += A(L1(j1),S1(k1),L2(j2),S2(k2))*B(L1(j1),L2(j2));
            }
        CHECK_CLOSE(res(S1(k1),S2(k2)),val,1E-5);
        }

    }

SECTION("NonContractProd")
    {

    IQTensor res = A / B;

    for(int j1 = 1; j1 <= L1.m(); ++j1)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
    for(int k1 = 1; k1 <= S1.m(); ++k1)
    for(int k2 = 1; k2 <= S2.m(); ++k2)
        {
        res(L1(j1),L2(j2),S1(k1),S2(k2));
        A(L1(j1),S1(k1),L2(j2),S2(k2));
        B(L1(j1),L2(j2));
        CHECK_CLOSE(res(L1(j1),L2(j2),S1(k1),S2(k2)),
                    A(L1(j1),S1(k1),L2(j2),S2(k2))*B(L1(j1),L2(j2)),1E-5);
        }

    }

SECTION("ComplexNonContractingProduct")
    {
    IQTensor Lr(L1(1),S1(2),L2(4)), Li(L1(1),S1(2),L2(4)),
            Rr(L1(1),S1(2),L2(4)), Ri(L1(1),S1(2),L2(4));

    Lr.randomize(); 
    Li.randomize(); 
    Rr.randomize();
    Ri.randomize();

    IQTensor L = Complex_1*Lr + Complex_i*Li;
    IQTensor R = Complex_1*Rr + Complex_i*Ri;

    IQTensor res1 = L / R;

    CHECK(hasindex(res1,L1));
    CHECK(hasindex(res1,S1));
    CHECK(hasindex(res1,L2));

    CHECK_EQUAL(res1.r(),3);

    ITensor resR(realPart(res1)),
            resI(imagPart(res1));

    ITensor rdiff = resR-(Lr/Rr-Li/Ri);
    ITensor idiff = resI-(Lr/Ri+Li/Rr);

    CHECK(rdiff.norm() < 1E-12);
    CHECK(idiff.norm() < 1E-12);
    }

SECTION("ITensorConversion")
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

///*
//SECTION("SymmetricDiag11")
//    {
//    IQTensor D,U;
//    IQIndex mid;
//    int mink,maxk;
//    C *= -2;
//    C.symmetricDiag11(L1,D,U,mid,mink,maxk);
//
//    IQTensor UD(U);
//    UD.prime(L1);
//    UD /= D;
//    ITensor diff = (dag(UD)*U - C).toITensor();
//    CHECK(diff.norm() < 1E-10);
//
//    IQTensor set1(dag(mid));
//    set1(mid(mink)) = 1;
//    U *= set1;
//    CHECK_CLOSE(D(mid(mink)),Dot(prime(U),C*U),1E-10);
//
//    }
//    */
//
SECTION("TieIndices")
    {
    IQTensor D1 = tieIndices(D,L1,prime(L1),L1);

    for(int k1 = 1; k1 <= L1.m(); ++k1)
    for(int k2 = 1; k2 <= L1.m(); ++k2)
    for(int k3 = 1; k3 <= S1.m(); ++k3)
        {
        CHECK_CLOSE(D1(L1(k1),prime(L1,2)(k2),S1(k3)),D(L1(k1),prime(L1)(k1),prime(L1,2)(k2),S1(k3)),1E-10);
        }
    }

SECTION("ToReal")
    {
    Real f = Global::random();
    IQTensor T(f);
    CHECK_CLOSE(T.toReal(),f,1E-5);

    //Default constructed IQTensor should
    //have a toReal value of zero
    IQTensor Z;
    CHECK_CLOSE(Z.toReal(),0,1E-5);
    }

SECTION("DotTest")
    {
    Real dotval1 = sqrt( Dot(dag(B),B) );
    //Dot should auto-fix arrows
    Real dotval2 = sqrt( Dot(B,B) );
    Real nval   = B.norm();
    CHECK_CLOSE(dotval1,nval,1E-5);
    CHECK_CLOSE(dotval2,nval,1E-5);
    }

SECTION("BraKetTest")
    {
    IQTensor R(L1(1),L2(1)),
             I(L1(1),L2(1));
    R.randomize();
    I.randomize();
    const Real rr = sqr(R.norm());
    const Real ii = sqr(I.norm());

    Complex z = BraKet(R,R);
    CHECK_CLOSE(z.real(),rr,1E-5);

    IQTensor T = Complex_1*R + Complex_i*I;
    z = BraKet(T,T);
    CHECK_CLOSE(z.real(),rr+ii,1E-5);
    CHECK(fabs(z.imag()) < 1E-12);

    z = BraKet(T,R);
    CHECK_CLOSE(z.real(),rr,1E-5);
    CHECK_CLOSE(z.imag(),-Dot(I,R),1E-5);

    z = BraKet(T,Complex_i*I);
    CHECK_CLOSE(z.real(),ii,1E-5);
    CHECK_CLOSE(z.imag(),Dot(I,R),1E-5);
    }

SECTION("Trace")
    {

    Real f = -Global::random();
    D *= f;

    IQTensor Dt = trace(D,dag(L1),prime(L1,2));

    for(int j2 = 1; j2 <= S1.m(); ++j2)
    for(int j1 = 1; j1 <= L1.m(); ++j1)
        {
        Real val = 0;
        for(int k1 = 1; k1 <= L1.m(); ++k1)
            {
            val += D(dag(L1)(k1),S1(j2),prime(L1)(j1),prime(L1,2)(k1));
            }
        CHECK_CLOSE(val,Dt(S1(j2),prime(L1)(j1)),1E-10);
        }

    IQTensor rho(L1(2),prime(L1)(2));
    rho.randomize();
    Real tr = trace(rho);
    rho /= tr;
    tr = trace(rho);
    CHECK(fabs(tr-1.) < 1E-11);
    }

SECTION("MapElems")
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

SECTION("Randomize")
    {
    IQTensor T(S1,S2,L2);
    T += ITensor(s1u,s2u,l2dd);
    T.randomize();
    //PrintDat(T);
    }

SECTION("RealImagPart")
    {
    IQTensor Z(dag(S1),prime(S1));
    Z(S1(1),prime(S1)(1)) = +1;
    Z(S1(2),prime(S1)(2)) = -1;

    IQTensor X(dag(S1),prime(S1));
    X(S1(1),prime(S1)(2)) = 1;
    X(S1(2),prime(S1)(1)) = 1;

    IQTensor ZiX = Complex_1*Z + Complex_i*X;
    IQTensor R(realPart(ZiX)),
             I(imagPart(ZiX));
    //PrintDat(R);
    //PrintDat(I);
    R -= Z;
    I -= X;
    CHECK_CLOSE(R.norm(),0,1E-5);
    CHECK_CLOSE(I.norm(),0,1E-5);

    //Test hc:

    ZiX.dag();
    R = realPart(ZiX);
    I = imagPart(ZiX);
    R -= Z;
    I += X;
    CHECK_CLOSE(R.norm(),0,1E-5);
    CHECK_CLOSE(I.norm(),0,1E-5);
    }

SECTION("ComplexMult")
    {
    IQTensor Z(dag(S1),prime(S1));
    Z(S1(1),prime(S1)(1)) = +1;
    Z(S1(2),prime(S1)(2)) = -1;

    IQTensor Y(dag(S1),prime(S1));
    Y(S1(1),prime(S1)(2)) =  1;
    Y(S1(2),prime(S1)(1)) = -1;
    Y *= Complex_i;

    //PrintDat(Y);

    IQTensor ZY = multSiteOps(Z,Y);
    //PrintDat(ZY);

    IQTensor YY = multSiteOps(Y,Y);
    //PrintDat(YY);
    }

SECTION("ComplexAdd")
    {
    IQTensor Z(dag(S1),prime(S1));
    Z(S1(1),prime(S1)(1)) = +1;
    Z(S1(2),prime(S1)(2)) = -1;

    IQTensor X(dag(S1),prime(S1));
    X(S1(1),prime(S1)(2)) = +1;
    X(S1(2),prime(S1)(1)) = +1;

    IQTensor iX = X * Complex_i;

    IQTensor R = Z + iX;

    CHECK_CLOSE((realPart(R)-Z).norm(),0,1E-5);
    CHECK_CLOSE((imagPart(R)-X).norm(),0,1E-5);

    }

SECTION("RandomizeTest")
    {
    IQTensor T(L1(1),S1(1),L2(4),S2(2));
    const QN D = div(T);
    T.randomize();
    CHECK_EQUAL(D,div(T));
    }

SECTION("Test_normLogNum")
    {
    IQTensor Z(dag(S1),prime(S1));
    ITensor blk1(s1u(1),prime(s1u)(1)),
            blk2(s1d(1),prime(s1d)(1));
    blk1 *= 0.1234;
    blk1 *= LogNumber(10,1);
    Z += blk1; 
    blk2 *= LogNumber(9,1);
    Z += blk2; 

    CHECK_CLOSE(Z.normLogNum().logNum(),log(sqrt(sqr(0.1234)*exp(20)+exp(18))),1E-5);

    }

SECTION("BigNorm")
    {
    IQTensor Z(dag(S1),prime(S1));
    ITensor blk1(s1u(1),prime(s1u)(1)),
            blk2(s1d(1),prime(s1d)(1));
    blk1 *= 0.1234;
    blk1 *= LogNumber(1000,1);
    Z += blk1; 
    blk2 *= LogNumber(999,1);
    Z += blk2; 

    //Mainly want to check that Z.normLogNum() doesn't overflow in this case
    CHECK_CLOSE(Z.normLogNum().logNum(),999.053,1E-3);
    }

SECTION("AddBlock")
    {
    ITensor b1(L1(1).indexqn(),L2(3).indexqn()),
            b2(L1(1).indexqn(),L2(2).indexqn());

    b1.randomize();
    b2.randomize();

    B += b1; //shouldn't throw

    CHECK_THROWS_AS(B += b2,ITError);
    
    }

SECTION("ComplexConvert")
    {
    IQTensor R(S1(1),L1(3)),
             I(S1(2),L1(1));
    R.randomize();
    I.randomize();
    R *= 0.1242;
    I *= -2.333;
    IQTensor T = R+Complex_i*I;

    //Global::debug1() = true;

    ITensor r = R.toITensor(),
            i = I.toITensor();
    ITensor t = T.toITensor();

    //Global::debug1() = false;


    CHECK((realPart(t)-r).norm() < 1E-12);
    CHECK((imagPart(t)-i).norm() < 1E-12);
    }

}
