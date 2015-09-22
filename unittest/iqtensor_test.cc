#include "test.h"
#include "itensor/iqtensor.h"
#include "itensor/util/set_scoped.h"
#include "itensor/util/count.h"

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

phi = randomTensor(S1(1),S2(1),L2(3));

A = randomTensor(L1(1),S1(1),L2(4),S2(2));

B = randomTensor(L1(1),L2(3));

C = randomTensor(dag(L1)(5),prime(L1)(5));

D = randomTensor(dag(L1)(3),S1(1),prime(L1)(3),prime(L1,2)(5));

SECTION("Boolean")
    {
    IQTensor t1;

    CHECK(!t1);

    CHECK(A);
    CHECK(B);
    CHECK(C);
    CHECK(D);
    }

SECTION("Contracting Product")
    {
    SECTION("Case 1")
        {
        auto R = A * dag(B);

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
            CHECK_CLOSE(R.real(S1(k1),S2(k2)),val);
            }
        }

    SECTION("Case 2")
        {
        auto R = A * dag(prime(B,L1));

        CHECK(!hasindex(R,L2));
        CHECK(hasindex(R,L1));
        CHECK(hasindex(R,S1));
        CHECK(hasindex(R,S2));

        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j1p = 1; j1p <= L1.m(); ++j1p)
            {
            Real val = 0;
            for(int j2 = 1; j2 <= L2.m(); ++j2)
                {
                val += A.real(L1(j1),S1(k1),L2(j2),S2(k2))*B.real(L1(j1p),L2(j2));
                }
            CHECK_CLOSE(R.real(prime(L1)(j1p),L1(j1),S1(k1),S2(k2)),val);
            }
        }

    SECTION("Regression Test 1")
        {
        auto s = IQIndex("S=1 site",
                  Index("Up",1,Site),QN(+2,0),
                  Index("Z0",1,Site),QN( 0,0),
                  Index("Dn",1,Site),QN(-2,0));

        auto sP = prime(s);

        IQIndexVal Up(s(1)),
                   UpP(sP(1)),
                   Dn(s(s.m())),
                   DnP(sP(s.m())),
                   Z0(s(2)),
                   Z0P(sP(2));

        IQTensor Op(dag(s),sP);
        Op.set(Z0,UpP,Sqrt2);
        Op.set(Dn,Z0P,Sqrt2);

        Op *= 0.5;
        //Op.scaleTo(1.); //This fixes the bug

        auto l0 = IQIndex("L0",Index("l0",1),QN());
        auto l1 = IQIndex("L1",Index("l1",3),QN());
        auto t = IQTensor(l0(1),l1(3));

        auto R = Op * t;

        for(auto i : count1(s.m()))
        for(auto iP : count1(sP.m()))
        for(auto j : count1(l1.m()))
            {
            auto val = Op.real(s(i),sP(iP)) * t.real(l0(1),l1(j));
            CHECK_CLOSE(val, R.real(s(i),sP(iP),l0(1),l1(j)) );
            }
        }

    }

SECTION("Addition and Subtraction")
    {
    SECTION("Case 1")
        {
        auto T1 = randomTensor(QN(0),L1,S1,L2,S2),
             T2 = randomTensor(QN(0),L1,S1,L2,S2);
        auto R = T1+T2;

        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
            {
            auto val = T1.real(L1(j1),S1(k1),L2(j2),S2(k2))+T2.real(L1(j1),S1(k1),L2(j2),S2(k2));
            CHECK_CLOSE(val,R.real(L1(j1),S1(k1),L2(j2),S2(k2)));
            }
        }

    SECTION("Case 2")
        {
        auto T1 = randomTensor(QN(0),L1,S1,L2,S2),
             T2 = randomTensor(QN(0),S1,S2,L1,L2);
        auto R = T1+T2;

        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
            {
            auto val = T1.real(L1(j1),S1(k1),L2(j2),S2(k2))+T2.real(L1(j1),S1(k1),L2(j2),S2(k2));
            CHECK_CLOSE(val,R.real(L1(j1),S1(k1),L2(j2),S2(k2)));
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
        CHECK_CLOSE( f( B.real(L1(j1),L2(j2)) ), 
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
            CHECK_CLOSE(phi.real(S1(k1),S2(k2),L2(j2)),itphi.real(Index(S1)(k1),Index(S2)(k2),Index(L2)(j2)));
        }

    SECTION("Case 2")
        {
        ITensor itA = A; //implicit conversion 
        for(int k1 = 1; k1 <= S1.m(); ++k1)
        for(int k2 = 1; k2 <= S2.m(); ++k2)
        for(int j1 = 1; j1 <= L1.m(); ++j1)
        for(int j2 = 1; j2 <= L2.m(); ++j2)
            CHECK_CLOSE(A.real(S1(k1),S2(k2),L1(j1),L2(j2)),
                        itA.real(Index(S1)(k1),Index(S2)(k2),Index(L1)(j1),Index(L2)(j2)));
        }
    }

SECTION("Combiner")
    {
    SECTION("Simple rank 2 combiner")
        {
        //Rank 2 combiner just replaces index
        auto s = IQIndex("s",
                         Index("s+1",2),QN(+1),
                         Index("s_0",2),QN( 0),
                         Index("s-1",2),QN(-1));
        auto C = combiner(s);

        auto T0 = randomTensor(QN(0),s,prime(s));
        auto R0 = C * T0;
        CHECK(norm(R0) == norm(T0));
        CHECK(div(R0) == div(T0));

        auto Tp1 = randomTensor(QN(+1),s,prime(s));
        auto Rp1 = C * Tp1;
        CHECK(norm(Rp1) == norm(Tp1));
        CHECK(div(Rp1) == div(Tp1));

        auto Tm1 = randomTensor(QN(-1),s,prime(s));
        auto Rm1 = C * Tm1;
        CHECK(norm(Rm1) == norm(Tm1));
        CHECK(div(Rm1) == div(Tm1));
        }
    SECTION("Combine / Uncombine 0 - No Permute")
        {
        auto T = randomTensor(QN(),L1,L2);
        auto C = combiner(L1);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        CHECK(ci); //check that ci was found
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
            {
            CHECK_CLOSE( T.real(L1(i1),L2(i2)), R.real(L1(i1),L2(i2)) );
            }
        }
    SECTION("Combine / Uncombine 1 - No Permute")
        {
        auto T = randomTensor(QN(),L1,L2);
        auto C = combiner(L1,L2);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        CHECK(ci);
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
            {
            CHECK_CLOSE( T.real(L1(i1),L2(i2)), R.real(L1(i1),L2(i2)) );
            }
        }

    SECTION("Combine / Uncombine 2 - No Permute")
        {
        auto T = randomTensor(QN(),L1,L2,S1);
        auto C = combiner(L1,L2);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        CHECK(hasindex(R,S1));
        CHECK(!hasindex(R,L1));
        CHECK(!hasindex(R,L2));
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        CHECK(!hasindex(R,ci));
        CHECK(hasindex(R,L1));
        CHECK(hasindex(R,L2));
        CHECK(hasindex(R,S1));
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
        for(int j1 = 1; j1 <= S1.m(); ++j1)
            {
            CHECK_CLOSE( T.real(L1(i1),L2(i2),S1(j1)), R.real(L1(i1),L2(i2),S1(j1)) );
            }
        }

    SECTION("Combine / Uncombine 3 - No Permute")
        {
        auto T = randomTensor(QN(),L1,S1,L2,S2);
        auto C = combiner(L1,S1);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
        for(int j1 = 1; j1 <= S1.m(); ++j1)
        for(int j2 = 1; j2 <= S2.m(); ++j2)
            {
            CHECK_CLOSE( T.real(L1(i1),L2(i2),S1(j1),S2(j2)), R.real(L1(i1),L2(i2),S1(j1),S2(j2)) );
            }
        }

    SECTION("Combine / Uncombine 4 - Permute")
        {
        auto T = randomTensor(QN(),L1,L2,S1,S2);
        auto C = combiner(L1,S1);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
        for(int j1 = 1; j1 <= S1.m(); ++j1)
        for(int j2 = 1; j2 <= S2.m(); ++j2)
            {
            CHECK_CLOSE( T.real(L1(i1),L2(i2),S1(j1),S2(j2)), R.real(L1(i1),L2(i2),S1(j1),S2(j2)) );
            }
        }

    SECTION("Combine / Uncombine 5 - Permute")
        {
        auto T = randomTensor(QN(),L1,L2,S1,S2);
        auto C = combiner(L1,S1,L2);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        CHECK(!hasindex(R,L1));
        CHECK(!hasindex(R,S1));
        CHECK(!hasindex(R,L2));
        CHECK(hasindex(R,S2));
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R = dag(C)*R; //uncombine
        CHECK(!hasindex(R,ci));
        CHECK(hasindex(R,L1));
        CHECK(hasindex(R,L2));
        CHECK(hasindex(R,S1));
        CHECK(hasindex(R,S2));
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
        for(int j1 = 1; j1 <= S1.m(); ++j1)
        for(int j2 = 1; j2 <= S2.m(); ++j2)
            {
            CHECK_CLOSE( T.real(L1(i1),L2(i2),S1(j1),S2(j2)), R.real(L1(i1),L2(i2),S1(j1),S2(j2)) );
            }
        }

    SECTION("Combine / Uncombine 6 - Permute")
        {
        auto T = randomTensor(QN(),L1,L2,S1,S2);
        auto C = combiner(S2,S1);
        auto R = C*T;
        auto ci = commonIndex(R,C); //get combined index
        //check that ci exists
        CHECK(ci);
        CHECK(hasindex(R,L1));
        CHECK(hasindex(R,L2));
        CHECK(!hasindex(R,S1));
        CHECK(!hasindex(R,S2));
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.m(); ++i1)
        for(int i2 = 1; i2 <= L2.m(); ++i2)
        for(int j1 = 1; j1 <= S1.m(); ++j1)
        for(int j2 = 1; j2 <= S2.m(); ++j2)
            {
            CHECK_CLOSE( T.real(L1(i1),L2(i2),S1(j1),S2(j2)), R.real(L1(i1),L2(i2),S1(j1),S2(j2)) );
            }
        }

    SECTION("Fragmented IQIndex Combiner Test 1")
        {
        auto i1 = IQIndex("i1",
                          Index("i1",2),QN(+1),
                          Index("i1",2),QN(+0),
                          Index("i1",2),QN(-1),
                          Index("i1",2),QN(+2),
                          Index("i1",2),QN(-2),
                          Index("i1",2),QN(+0),
                          Index("i1",2),QN(-1));

        auto i2 = IQIndex("i2",
                          Index("i2",2),QN(+1),
                          Index("i2",4),QN(-1),
                          Index("i2",2),QN(+0),
                          Index("i2",3),QN(+2),
                          Index("i2",2),QN(+0),
                          Index("i2",2),QN(-1),
                          Index("i2",3),QN(+1));

        auto i3 = IQIndex("i3",
                          Index("i3",3),QN(-1),
                          Index("i3",2),QN(+0),
                          Index("i3",2),QN(+1),
                          Index("i3",2),QN(+0),
                          Index("i3",4),QN(-1),
                          Index("i3",2),QN(+0),
                          Index("i3",2),QN(+1));

        auto flux = QN(-2);
        auto T = randomTensor(flux,i1,prime(i3),i2,i3,prime(i2));
        auto C = combiner(i1,i2,prime(i2));
        auto R = C * T;
        auto nT = dag(C) * R;

        CHECK(div(T) == div(R));
        CHECK(div(T) == div(nT));
        CHECK(norm(T-nT) < 1E-11);
        }

    SECTION("Fragmented IQIndex Combiner Test 2")
        {
        auto i1 = IQIndex("i1",
                          Index("i1",2),QN(+1),
                          Index("i1",2),QN(+0),
                          Index("i1",2),QN(-1),
                          Index("i1",2),QN(+2),
                          Index("i1",2),QN(-2),
                          Index("i1",2),QN(+0),
                          Index("i1",2),QN(-1));

        auto i2 = IQIndex("i2",
                          Index("i2",2),QN(+1),
                          Index("i2",4),QN(-1),
                          Index("i2",2),QN(+0),
                          Index("i2",3),QN(+2),
                          Index("i2",2),QN(+0),
                          Index("i2",2),QN(-1),
                          Index("i2",3),QN(+1));

        auto i3 = IQIndex("i3",
                          Index("i3",3),QN(-1),
                          Index("i3",2),QN(+0),
                          Index("i3",2),QN(+1),
                          Index("i3",2),QN(+0),
                          Index("i3",4),QN(-1),
                          Index("i3",2),QN(+0),
                          Index("i3",2),QN(+1));

        auto flux = QN(-2);
        auto T = randomTensor(flux,i1,prime(i3),i2,i3,prime(i2));
        auto C = combiner(i3,i1);
        auto R = C * T;
        auto nT = dag(C) * R;

        CHECK(div(T) == div(R));
        CHECK(div(T) == div(nT));
        CHECK(norm(T-nT) < 1E-11);
        }

    }

SECTION("Scalar")
    {
    auto T1 = randomTensor(QN(),L1,L2,S1,S2);
    auto T2 = randomTensor(QN(),S1,L2,S2,L1);
    auto S = T1*dag(T2);
    CHECK(S.r() == 0);

    Real val = 0;
    for(int i1 = 1; i1 <= L1.m(); ++i1)
    for(int i2 = 1; i2 <= L2.m(); ++i2)
    for(int j1 = 1; j1 <= S1.m(); ++j1)
    for(int j2 = 1; j2 <= S2.m(); ++j2)
        {
        val += T1.real(L1(i1),L2(i2),S1(j1),S2(j2))*T2.real(L1(i1),L2(i2),S1(j1),S2(j2));
        }

    CHECK_CLOSE(val, S.real());
    CHECK_CLOSE(fabs(val), norm(S));
    }



//SECTION("TieIndices")
//    {
//    IQTensor D1 = tieIndices(D,L1,prime(L1),L1);
//
//    for(int k1 = 1; k1 <= L1.m(); ++k1)
//    for(int k2 = 1; k2 <= L1.m(); ++k2)
//    for(int k3 = 1; k3 <= S1.m(); ++k3)
//        {
//        CHECK_DIFF(D1(L1(k1),prime(L1,2)(k2),S1(k3)),D(L1(k1),prime(L1)(k1),prime(L1,2)(k2),S1(k3)),1E-10);
//        }
//    }


//SECTION("DotTest")
//    {
//    Real dotval1 = sqrt( Dot(dag(B),B) );
//    //Dot should auto-fix arrows
//    Real dotval2 = sqrt( Dot(B,B) );
//    Real nval   = B.norm();
//    CHECK_DIFF(dotval1,nval,1E-5);
//    CHECK_DIFF(dotval2,nval,1E-5);
//    }
//
//SECTION("BraKetTest")
//    {
//    auto R = randomTensor(L1(1),L2(1)),
//         I = randomTensor(L1(1),L2(1));
//    const Real rr = sqr(R.norm());
//    const Real ii = sqr(I.norm());
//
//    Complex z = BraKet(R,R);
//    CHECK_DIFF(z.real(),rr,1E-5);
//
//    IQTensor T = Complex_1*R + Complex_i*I;
//    z = BraKet(T,T);
//    CHECK_DIFF(z.real(),rr+ii,1E-5);
//    CHECK(fabs(z.imag()) < 1E-12);
//
//    z = BraKet(T,R);
//    CHECK_DIFF(z.real(),rr,1E-5);
//    CHECK_DIFF(z.imag(),-Dot(I,R),1E-5);
//
//    z = BraKet(T,Complex_i*I);
//    CHECK_DIFF(z.real(),ii,1E-5);
//    CHECK_DIFF(z.imag(),Dot(I,R),1E-5);
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
//        CHECK_DIFF(val,Dt(S1(j2),prime(L1)(j1)),1E-10);
//        }
//
//    auto rho = randomTensor(L1(2),prime(L1)(2));
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
//    CHECK_DIFF(R.norm(),0,1E-5);
//    CHECK_DIFF(I.norm(),0,1E-5);
//
//    //Test hc:
//
//    ZiX.dag();
//    R = realPart(ZiX);
//    I = imagPart(ZiX);
//    R -= Z;
//    I += X;
//    CHECK_DIFF(R.norm(),0,1E-5);
//    CHECK_DIFF(I.norm(),0,1E-5);
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
//    CHECK_DIFF((realPart(R)-Z).norm(),0,1E-5);
//    CHECK_DIFF((imagPart(R)-X).norm(),0,1E-5);
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
//    CHECK_DIFF(Z.normLogNum().logNum(),log(sqrt(sqr(0.1234)*exp(20)+exp(18))),1E-5);
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
//    CHECK_DIFF(Z.normLogNum().logNum(),999.053,1E-3);
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
//    auto R = randomTensor(S1(1),L1(3)),
//         I = randomTensor(S1(2),L1(1));
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
