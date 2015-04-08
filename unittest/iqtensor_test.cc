#include "test.h"
#include "iqtensor.h"

using namespace itensor;
using namespace std;

struct FuncObj
    {
    template<typename T>
    T
    operator()(T x) const { return x*x; }
    };


TEST_CASE("IQTensorTest")
{

Index s1u("S1+1",1,Site);
Index s1d("S1-1",1,Site);
Index s2u("S2+1",1,Site);
Index s2d("S2-1",1,Site);
Index l1u("L1+1",2,Link);
Index l10("L1_0",2,Link);
Index l1d("L1-1",2,Link);
Index l2uu("L2+2",2,Link);
Index l20("L2_0",2,Link);
Index l2dd("L2-2",2,Link);

IQIndex S1,S2,L1,L2;

IQTensor phi,A,B,C,D;

S1 = IQIndex("S1",Out,
             s1u,QN(+1),
             s1d,QN(-1));
S2 = IQIndex("S2",Out,
             s2u,QN(+1),
             s2d,QN(-1));
L1 = IQIndex("L1",Out,
             l1u,QN(+1),
             l10,QN( 0),
             l1d,QN(-1));
             
L2 = IQIndex("L2",Out,
             l2uu,QN(+2),
             l20,QN( 0),
             l2dd,QN(-2));

phi = randIQT(S1(1),S2(1),L2(3));

A = randIQT(L1(1),S1(1),L2(4),S2(2));

B = randIQT(L1(1),L2(3));

C = randIQT(dag(L1)(5),prime(L1)(5));

D = randIQT(dag(L1)(3),S1(1),prime(L1)(3),prime(L1,2)(5));

SECTION("Boolean")
    {
    IQTensor t1;

    CHECK(!t1);

    CHECK(A);
    CHECK(B);
    CHECK(C);
    CHECK(D);
    }

//SECTION("Constructors")
//    {
//    Real f = Global::random();
//    IQTensor rZ(f);
//
//    CHECK_EQUAL(rZ.r(),0);
//    CHECK_REQUAL(norm(rZ),f);
//    }

//SECTION("ToReal")
//    {
//    Real f = Global::random();
//    IQTensor T(f);
//    PrintData(T);
//    CHECK_REQUAL(T.real(),f);
//    }

SECTION("Contracting Product")
    {
    SECTION("Case 1")
        {
        //PrintData(A);
        //PrintData(B);
        auto R = A*dag(B);
        //PrintData(R);

        CHECK(hasindex(R,S1));
        CHECK(hasindex(R,S2));
        CHECK(!hasindex(R,L1));
        CHECK(!hasindex(R,L2));

        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
            {
            Real val = 0;
            for(int j1 = 1; j1 <= L1.m(); ++j1)
            for(int j2 = 1; j2 <= L2.m(); ++j2)
                {
                //printf("val += %f*%f",A.real(L1(j1),S1(k1),L2(j2),S2(k2)),B.real(L1(j1),L2(j2)));
                val += A.real(L1(j1),S1(k1),L2(j2),S2(k2))*B.real(L1(j1),L2(j2));
                //printfln(" (now val=%f)",val);
                }
            //printfln("val = %f, R.real(S1(%d),S2(%d))=%f",val,k1,k2,R.real(S1(k1),S2(k2)));
            CHECK_REQUAL(R.real(S1(k1),S2(k2)),val);
            }
        }

    }

SECTION("Addition and Subtraction")
    {
    SECTION("Case 1")
        {
        auto T1 = randIQT(QN(0),L1,S1,L2,S2),
             T2 = randIQT(QN(0),L1,S1,L2,S2);
        auto R = T1+T2;

        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
            {
            auto val = T1.real(L1(j1),S1(k1),L2(j2),S2(k2))+T2.real(L1(j1),S1(k1),L2(j2),S2(k2));
            CHECK_REQUAL(val,R.real(L1(j1),S1(k1),L2(j2),S2(k2)));
            }
        }

    SECTION("Case 2")
        {
        auto T1 = randIQT(QN(0),L1,S1,L2,S2),
             T2 = randIQT(QN(0),S1,S2,L1,L2);
        auto R = T1+T2;

        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
            {
            auto val = T1.real(L1(j1),S1(k1),L2(j2),S2(k2))+T2.real(L1(j1),S1(k1),L2(j2),S2(k2));
            CHECK_REQUAL(val,R.real(L1(j1),S1(k1),L2(j2),S2(k2)));
            }
        }

    }

SECTION("Apply")
    {
    IQTensor B1(B);

    FuncObj f;
    B1.apply(f);

    for(int j1 = 1; j1 <= L1.m(); ++j1)
    for(int j2 = 1; j2 <= L2.m(); ++j2)
        {
        CHECK_REQUAL( f( B.real(L1(j1),L2(j2)) ), 
                    B1.real(L1(j1),L2(j2)));
        }
    }

SECTION("RandomizeTest")
    {
    IQTensor T(L1(1),S1(1),L2(4),S2(2));
    const QN D = div(T);
    T = randomize(T);
    CHECK_EQUAL(D,div(T));
    }

SECTION("ITensor Conversion")
    {
    SECTION("Case 1")
        {
        auto itphi = toITensor(phi);
        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
            CHECK_REQUAL(phi.real(S1(k1),S2(k2),L2(j2)),itphi.real(Index(S1)(k1),Index(S2)(k2),Index(L2)(j2)));
        }

    SECTION("Case 2")
        {
        ITensor itA = A; //implicit conversion 
        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
            CHECK_REQUAL(A.real(S1(k1),S2(k2),L1(j1),L2(j2)),
                        itA.real(Index(S1)(k1),Index(S2)(k2),Index(L1)(j1),Index(L2)(j2)));
        }
    }

SECTION("Combiner")
    {
    SECTION("Combine / Uncombine 0 - No Permute")
        {
        auto T = randIQT(QN(),L1,L2);
        auto C = combiner(L1);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        //check that all elements of T accounted for in R
        CHECK_REQUAL(norm(T),norm(R));
        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
            {
            CHECK_REQUAL( T.real(L1(i1),L2(i2)), R.real(L1(i1),L2(i2)) );
            }
        }
    SECTION("Combine / Uncombine 1 - No Permute")
        {
        auto T = randIQT(QN(),L1,L2);
        auto C = combiner(L1,L2);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        //check that all elements of T accounted for in R
        CHECK_REQUAL(norm(T),norm(R));
        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
            {
            CHECK_REQUAL( T.real(L1(i1),L2(i2)), R.real(L1(i1),L2(i2)) );
            }
        }

    SECTION("Combine / Uncombine 2 - No Permute")
        {
        auto T = randIQT(QN(),L1,L2,S1);
        auto C = combiner(L1,L2);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        //check that all elements of T accounted for in R
        CHECK_REQUAL(norm(T),norm(R));
        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
        for(int j1 = 1; j1 <= S1.m(); ++j1)
            {
            CHECK_REQUAL( T.real(L1(i1),L2(i2),S1(j1)), R.real(L1(i1),L2(i2),S1(j1)) );
            }
        }

    SECTION("Combine / Uncombine 3 - No Permute")
        {
        auto T = randIQT(QN(),L1,S1,L2,S2);
        auto C = combiner(L1,S1);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        //check that all elements of T accounted for in R
        CHECK_REQUAL(norm(T),norm(R));
        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
        for(int j1 = 1; j1 <= S1.m(); ++j1)
        for(int j2 = 1; j2 <= S2.m(); ++j2)
            {
            CHECK_REQUAL( T.real(L1(i1),L2(i2),S1(j1),S2(j2)), R.real(L1(i1),L2(i2),S1(j1),S2(j2)) );
            }
        }

    }



//SECTION("TieIndices")
//    {
//    IQTensor D1 = tieIndices(D,L1,prime(L1),L1);
//
//    for(int k1 = 1; k1 <= L1.m(); ++k1)
//    for(int k2 = 1; k2 <= L1.m(); ++k2)
//    for(int k3 = 1; k3 <= S1.m(); ++k3)
//        {
//        CHECK_CLOSE(D1(L1(k1),prime(L1,2)(k2),S1(k3)),D(L1(k1),prime(L1)(k1),prime(L1,2)(k2),S1(k3)),1E-10);
//        }
//    }


//SECTION("DotTest")
//    {
//    Real dotval1 = sqrt( Dot(dag(B),B) );
//    //Dot should auto-fix arrows
//    Real dotval2 = sqrt( Dot(B,B) );
//    Real nval   = B.norm();
//    CHECK_CLOSE(dotval1,nval,1E-5);
//    CHECK_CLOSE(dotval2,nval,1E-5);
//    }
//
//SECTION("BraKetTest")
//    {
//    auto R = randIQT(L1(1),L2(1)),
//         I = randIQT(L1(1),L2(1));
//    const Real rr = sqr(R.norm());
//    const Real ii = sqr(I.norm());
//
//    Complex z = BraKet(R,R);
//    CHECK_CLOSE(z.real(),rr,1E-5);
//
//    IQTensor T = Complex_1*R + Complex_i*I;
//    z = BraKet(T,T);
//    CHECK_CLOSE(z.real(),rr+ii,1E-5);
//    CHECK(fabs(z.imag()) < 1E-12);
//
//    z = BraKet(T,R);
//    CHECK_CLOSE(z.real(),rr,1E-5);
//    CHECK_CLOSE(z.imag(),-Dot(I,R),1E-5);
//
//    z = BraKet(T,Complex_i*I);
//    CHECK_CLOSE(z.real(),ii,1E-5);
//    CHECK_CLOSE(z.imag(),Dot(I,R),1E-5);
//    }

//SECTION("Trace")
//    {
//
//    Real f = -Global::random();
//    D *= f;
//
//    IQTensor Dt = trace(D,dag(L1),prime(L1,2));
//
//    for(int j2 = 1; j2 <= S1.m(); ++j2)
//    for(int j1 = 1; j1 <= L1.m(); ++j1)
//        {
//        Real val = 0;
//        for(int k1 = 1; k1 <= L1.m(); ++k1)
//            {
//            val += D(dag(L1)(k1),S1(j2),prime(L1)(j1),prime(L1,2)(k1));
//            }
//        CHECK_CLOSE(val,Dt(S1(j2),prime(L1)(j1)),1E-10);
//        }
//
//    auto rho = randIQT(L1(2),prime(L1)(2));
//    Real tr = trace(rho);
//    rho /= tr;
//    tr = trace(rho);
//    CHECK(fabs(tr-1.) < 1E-11);
//    }


//SECTION("RealImagPart")
//    {
//    IQTensor Z(dag(S1),prime(S1));
//    Z(S1(1),prime(S1)(1)) = +1;
//    Z(S1(2),prime(S1)(2)) = -1;
//
//    IQTensor X(dag(S1),prime(S1));
//    X(S1(1),prime(S1)(2)) = 1;
//    X(S1(2),prime(S1)(1)) = 1;
//
//    IQTensor ZiX = Complex_1*Z + Complex_i*X;
//    IQTensor R(realPart(ZiX)),
//             I(imagPart(ZiX));
//    //PrintDat(R);
//    //PrintDat(I);
//    R -= Z;
//    I -= X;
//    CHECK_CLOSE(R.norm(),0,1E-5);
//    CHECK_CLOSE(I.norm(),0,1E-5);
//
//    //Test hc:
//
//    ZiX.dag();
//    R = realPart(ZiX);
//    I = imagPart(ZiX);
//    R -= Z;
//    I += X;
//    CHECK_CLOSE(R.norm(),0,1E-5);
//    CHECK_CLOSE(I.norm(),0,1E-5);
//    }

//SECTION("ComplexMult")
//    {
//    IQTensor Z(dag(S1),prime(S1));
//    Z(S1(1),prime(S1)(1)) = +1;
//    Z(S1(2),prime(S1)(2)) = -1;
//
//    IQTensor Y(dag(S1),prime(S1));
//    Y(S1(1),prime(S1)(2)) =  1;
//    Y(S1(2),prime(S1)(1)) = -1;
//    Y *= Complex_i;
//
//    //PrintDat(Y);
//
//    IQTensor ZY = multSiteOps(Z,Y);
//    //PrintDat(ZY);
//
//    IQTensor YY = multSiteOps(Y,Y);
//    //PrintDat(YY);
//    }

//SECTION("ComplexAdd")
//    {
//    IQTensor Z(dag(S1),prime(S1));
//    Z(S1(1),prime(S1)(1)) = +1;
//    Z(S1(2),prime(S1)(2)) = -1;
//
//    IQTensor X(dag(S1),prime(S1));
//    X(S1(1),prime(S1)(2)) = +1;
//    X(S1(2),prime(S1)(1)) = +1;
//
//    IQTensor iX = X * Complex_i;
//
//    IQTensor R = Z + iX;
//
//    CHECK_CLOSE((realPart(R)-Z).norm(),0,1E-5);
//    CHECK_CLOSE((imagPart(R)-X).norm(),0,1E-5);
//
//    }

//SECTION("Test_normLogNum")
//    {
//    IQTensor Z(dag(S1),prime(S1));
//    ITensor blk1(s1u(1),prime(s1u)(1)),
//            blk2(s1d(1),prime(s1d)(1));
//    blk1 *= 0.1234;
//    blk1 *= LogNumber(10,1);
//    Z += blk1; 
//    blk2 *= LogNumber(9,1);
//    Z += blk2; 
//
//    CHECK_CLOSE(Z.normLogNum().logNum(),log(sqrt(sqr(0.1234)*exp(20)+exp(18))),1E-5);
//
//    }

//SECTION("BigNorm")
//    {
//    IQTensor Z(dag(S1),prime(S1));
//    ITensor blk1(s1u(1),prime(s1u)(1)),
//            blk2(s1d(1),prime(s1d)(1));
//    blk1 *= 0.1234;
//    blk1 *= LogNumber(1000,1);
//    Z += blk1; 
//    blk2 *= LogNumber(999,1);
//    Z += blk2; 
//
//    //Mainly want to check that Z.normLogNum() doesn't overflow in this case
//    CHECK_CLOSE(Z.normLogNum().logNum(),999.053,1E-3);
//    }

//SECTION("AddBlock")
//    {
//    ITensor b1(L1(1).indexqn(),L2(3).indexqn()),
//            b2(L1(1).indexqn(),L2(2).indexqn());
//
//    b1.randomize();
//    b2.randomize();
//
//    B += b1; //shouldn't throw
//
//    CHECK_THROWS_AS(B += b2,ITError);
//    
//    }

//SECTION("ComplexConvert")
//    {
//    auto R = randIQT(S1(1),L1(3)),
//         I = randIQT(S1(2),L1(1));
//    R *= 0.1242;
//    I *= -2.333;
//    IQTensor T = R+Complex_i*I;
//
//    //Global::debug1() = true;
//
//    ITensor r = R.toITensor(),
//            i = I.toITensor();
//    ITensor t = T.toITensor();
//
//    //Global::debug1() = false;
//
//
//    CHECK((realPart(t)-r).norm() < 1E-12);
//    CHECK((imagPart(t)-i).norm() < 1E-12);
//    }

}
