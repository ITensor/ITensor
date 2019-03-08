#include "test.h"
#include "itensor/iqtensor.h"
#include "itensor/util/set_scoped.h"
#include "itensor/util/range.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using namespace std;

struct FuncObj
    {
    template<typename T>
    T
    operator()(T x) const { return x*x; }
    };

bool
isIQTensor(ITensor const& T) { return false; }
bool
isIQTensor(IQTensor const& T) { return true; }

enum class QType { 
            OtherType, 
            QDenseReal, 
            QDenseCplx, 
            QDiagReal, 
            QDiagRealAllSame, 
            QDiagCplx, 
            QDiagCplxAllSame, 
            QCombiner,
            DiagReal,
            DiagCplx
          };

struct GetQType
    {
    QType operator()(QDenseReal const& d) const { return QType::QDenseReal; }
    QType operator()(QDenseCplx const& d) const { return QType::QDenseCplx; }
    QType operator()(QDiagReal const& d) const { return d.allSame() ? QType::QDiagRealAllSame : QType::QDiagReal; }
    QType operator()(QDiagCplx const& d) const { return d.allSame() ? QType::QDiagCplxAllSame : QType::QDiagCplx; }
    QType operator()(QCombiner const& d) const { return QType::QCombiner; }
    QType operator()(DiagReal const& d) const { return QType::DiagReal; }
    QType operator()(DiagCplx const& d) const { return QType::DiagCplx; }
    template<typename T>
    QType operator()(T const& d) const { return QType::OtherType; }
    };

template<typename I>
QType
typeOf(ITensorT<I> const& t) 
    { 
    return applyFunc(GetQType{},t.store()); 
    }

std::ostream&
operator<<(std::ostream& s, QType t)
    {
    if(t == QType::OtherType) s << "OtherType";
    else if(t == QType::QDenseReal) s << "QDenseReal";
    else if(t == QType::QDenseCplx) s << "QDenseCplx";
    else if(t == QType::QDiagReal) s << "QDiagReal";
    else if(t == QType::QDiagRealAllSame) s << "QDiagRealAllSame";
    else if(t == QType::QDiagCplx) s << "QDiagCplx";
    else if(t == QType::QDiagCplxAllSame) s << "QDiagCplxAllSame";
    else if(t == QType::QCombiner) s << "QCombiner";
    else if(t == QType::DiagReal) s << "DiagReal";
    else if(t == QType::DiagCplx) s << "DiagCplx";
    else Error("Unrecognized QType value");
    return s;
    }

TEST_CASE("IQTensorTest")
{

auto S1 = IQIndex(Index(1,"S1,Site"),QN(-1),
                  Index(1,"S1,Site"),QN(+1));
auto S2 = IQIndex(Index(1,"S2,Site"),QN(-1),
                  Index(1,"S2,Site"),QN(+1));
auto S3 = IQIndex(Index(1,"S3,Site"),QN(-1),
                  Index(1,"S3,Site"),QN(+1));

auto S4 = IQIndex(Index(1,"S4,Site"),QN(-1),
                  Index(1,"S4,Site"),QN(+1));

auto L1 = IQIndex(Index(2,"L1,Link"),QN(+1),
                  Index(2,"L1,Link"),QN( 0),
                  Index(2,"L1,Link"),QN(-1));
             
auto L2 = IQIndex(Index(2,"L2,Link"),QN(+2),
                  Index(2,"L2,Link"),QN( 0),
                  Index(2,"L2,Link"),QN(-2));

auto phi = randomTensor(S1(1),S2(1),L2(3));
auto A = randomTensor(L1(1),S1(1),L2(4),S2(2));
auto B = randomTensor(L1(1),L2(3));
auto C = randomTensor(dag(L1)(5),prime(L1)(5));
auto D = randomTensor(dag(L1)(3),S1(1),prime(L1)(3),prime(L1,2)(5));

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

        CHECK(hasIndex(R,S1));
        CHECK(hasIndex(R,S2));
        CHECK(!hasIndex(R,L1));
        CHECK(!hasIndex(R,L2));

        for(int k1 = 1; k1 <= S1.dim(); ++k1)
        for(int k2 = 1; k2 <= S2.dim(); ++k2)
            {
            Real val = 0;
            for(int j1 = 1; j1 <= L1.dim(); ++j1)
            for(int j2 = 1; j2 <= L2.dim(); ++j2)
                {
                //printf("val += %f*%f",elt(A,L1(j1),S1(k1),L2(j2),S2(k2)),elt(B,L1(j1),L2(j2)));
                val += elt(A,L1(j1),S1(k1),L2(j2),S2(k2))*elt(B,L1(j1),L2(j2));
                //printfln(" (now val=%f)",val);
                }
            //printfln("val = %f, R.elt(S1(%d),S2(%d))=%f",val,k1,k2,R.elt(S1(k1),S2(k2)));
            CHECK_CLOSE(R.elt(S1(k1),S2(k2)),val);
            }
        }

    SECTION("Case 2")
        {
        auto R = A * dag(prime(B,L1));

        CHECK(!hasIndex(R,L2));
        CHECK(hasIndex(R,L1));
        CHECK(hasIndex(R,S1));
        CHECK(hasIndex(R,S2));

        for(int k1 = 1; k1 <= S1.dim(); ++k1)
        for(int k2 = 1; k2 <= S2.dim(); ++k2)
        for(int j1 = 1; j1 <= L1.dim(); ++j1)
        for(int j1p = 1; j1p <= L1.dim(); ++j1p)
            {
            Real val = 0;
            for(int j2 = 1; j2 <= L2.dim(); ++j2)
                {
                val += A.elt(L1(j1),S1(k1),L2(j2),S2(k2))*B.elt(L1(j1p),L2(j2));
                }
            CHECK_CLOSE(R.elt(prime(L1)(j1p),L1(j1),S1(k1),S2(k2)),val);
            }
        }

    SECTION("Regression Test 1")
        {
        auto s = IQIndex(Index(1,"S=1,Site"),QN("Sz=",+2),
                         Index(1,"S=1,Site"),QN("Sz=", 0),
                         Index(1,"S=1,Site"),QN("Sz=",-2));

        auto sP = prime(s);

        IQIndexVal Up(s(1)),
                   UpP(sP(1)),
                   Dn(s(s.dim())),
                   DnP(sP(s.dim())),
                   Z0(s(2)),
                   Z0P(sP(2));

        IQTensor Op(dag(s),sP);
        Op.set(Z0,UpP,Sqrt2);
        Op.set(Dn,Z0P,Sqrt2);

        Op *= 0.5;
        //Op.scaleTo(1.); //This fixes the bug

        auto l0 = IQIndex(Index(1,"L0"),QN());
        auto l1 = IQIndex(Index(3,"L1"),QN());
        auto t = setElt(l0(1),l1(3));

        auto R = Op * t;

        for(auto i : range1(s.dim()))
        for(auto iP : range1(sP.dim()))
        for(auto j : range1(l1.dim()))
            {
            auto val = Op.elt(s(i),sP(iP)) * t.elt(l0(1),l1(j));
            CHECK_CLOSE(val, R.elt(s(i),sP(iP),l0(1),l1(j)) );
            }
        }

    SECTION("Regression Test 2")
        {
        //Feb 10 2016: was getting an error when trying
		//to print the result C of the following contraction
        //Bug was that doTask(CalcDiv) was being too stingy about
        //QDense storage with no blocks and throwing an exception
		auto s = IQIndex(Index(1,"s,Site"),QN(1),Index(1,"s,Site"),QN(-1));
		auto l = IQIndex(Index(1,"l,Link"),QN(3));

		auto A = IQTensor(s,l);
		A.set(s(1),l(1),1);

		auto B = IQTensor(s,prime(l));
		B.set(s(2),prime(l)(1),1);

		auto C = A*dag(B);

		auto q = div(C);
		CHECK(q == QN());
        }

    }

SECTION("Addition and Subtraction")
    {
    SECTION("Case 1")
        {
        auto T1 = randomTensor(QN(0),L1,S1,L2,S2),
             T2 = randomTensor(QN(0),L1,S1,L2,S2);
        auto R = T1+T2;

        for(int j1 = 1; j1 <= L1.dim(); ++j1)
        for(int j2 = 1; j2 <= L2.dim(); ++j2)
        for(int k1 = 1; k1 <= S1.dim(); ++k1)
        for(int k2 = 1; k2 <= S2.dim(); ++k2)
            {
            auto val = T1.elt(L1(j1),S1(k1),L2(j2),S2(k2))+T2.elt(L1(j1),S1(k1),L2(j2),S2(k2));
            CHECK_CLOSE(val,R.elt(L1(j1),S1(k1),L2(j2),S2(k2)));
            }
        }

    SECTION("Case 2")
        {
        auto T1 = randomTensor(QN(0),L1,S1,L2,S2),
             T2 = randomTensor(QN(0),S1,S2,L1,L2);
        auto R = T1+T2;

        for(int j1 = 1; j1 <= L1.dim(); ++j1)
        for(int j2 = 1; j2 <= L2.dim(); ++j2)
        for(int k1 = 1; k1 <= S1.dim(); ++k1)
        for(int k2 = 1; k2 <= S2.dim(); ++k2)
            {
            auto val = T1.elt(L1(j1),S1(k1),L2(j2),S2(k2))+T2.elt(L1(j1),S1(k1),L2(j2),S2(k2));
            CHECK_CLOSE(val,R.elt(L1(j1),S1(k1),L2(j2),S2(k2)));
            }
        }

    }

SECTION("Apply")
    {
    IQTensor B1(B);

    FuncObj f;
    B1.apply(f);

    for(int j1 = 1; j1 <= L1.dim(); ++j1)
    for(int j2 = 1; j2 <= L2.dim(); ++j2)
        {
        CHECK_CLOSE( f( B.elt(L1(j1),L2(j2)) ), 
                    B1.elt(L1(j1),L2(j2)));
        }
    }

SECTION("RandomizeTest")
    {
    auto T = setElt(L1(1),S1(1),L2(4),S2(2));
    const QN D = div(T);
    T.randomize();
    CHECK_EQUAL(D,div(T));
    }

SECTION("QDense ITensor Conversion")
    {
    SECTION("Case 1")
        {
        auto itphi = toITensor(phi);
        for(int k1 = 1; k1 <= S1.dim(); ++k1)
        for(int k2 = 1; k2 <= S2.dim(); ++k2)
        for(int j2 = 1; j2 <= L2.dim(); ++j2)
            CHECK_CLOSE(phi.elt(S1(k1),S2(k2),L2(j2)),itphi.elt(Index(S1)(k1),Index(S2)(k2),Index(L2)(j2)));
        }

    SECTION("Case 2")
        {
        ITensor itA = A; //implicit conversion 
        for(int k1 = 1; k1 <= S1.dim(); ++k1)
        for(int k2 = 1; k2 <= S2.dim(); ++k2)
        for(int j1 = 1; j1 <= L1.dim(); ++j1)
        for(int j2 = 1; j2 <= L2.dim(); ++j2)
            CHECK_CLOSE(A.elt(S1(k1),S2(k2),L1(j1),L2(j2)),
                        itA.elt(Index(S1)(k1),Index(S2)(k2),Index(L1)(j1),Index(L2)(j2)));
        }
    }

SECTION("QDiag ITensor Conversion")
    {
    SECTION("Case 1 (delta/allSame)")
        {
        auto D = delta(dag(L1),prime(L1));
        auto d = toITensor(D);
        CHECK(typeOf(d) == QType::DiagReal);
        Index l1 = L1;
        for(auto i : range1(L1))
        for(auto j : range1(L1))
            {
            CHECK_CLOSE(D.elt(L1(i),prime(L1)(j)),d.elt(l1(i),prime(l1)(j)));
            }
        }
    SECTION("Case 2")
        {
        auto D = delta(dag(L1),prime(L1));
        D.randomize();
        auto d = toITensor(D);
        CHECK(typeOf(d) == QType::DiagReal);
        Index l1 = L1;
        for(auto i : range1(L1))
        for(auto j : range1(L1))
            {
            CHECK_CLOSE(D.elt(L1(i),prime(L1)(j)),d.elt(l1(i),prime(l1)(j)));
            }
        }
    }

SECTION("Combiner")
    {
    SECTION("Simple order 2 combiner")
        {
        //Rank 2 combiner just replaces index
        auto s = IQIndex(Index(2),QN(+1),
                         Index(2),QN( 0),
                         Index(2),QN(-1));
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
        for(int i1 = 1; i1 <= L1.dim(); ++i1)
        for(int i2 = 1; i2 <= L2.dim(); ++i2)
            {
            CHECK_CLOSE( elt(T,L1(i1),L2(i2)), elt(R,L1(i1),L2(i2)) );
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
        for(int i1 = 1; i1 <= L1.dim(); ++i1)
        for(int i2 = 1; i2 <= L2.dim(); ++i2)
            {
            CHECK_CLOSE( T.elt(L1(i1),L2(i2)), R.elt(L1(i1),L2(i2)) );
            }
        }

    SECTION("Combine / Uncombine 2 - No Permute")
        {
        auto T = randomTensor(QN(),L1,L2,S1);
        auto C = combiner(L1,L2);
        auto R = T*C;
        auto ci = commonIndex(R,C); //get combined index
        CHECK(ci); //check that ci exists
        CHECK(hasIndex(R,S1));
        CHECK(!hasIndex(R,L1));
        CHECK(!hasIndex(R,L2));
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        CHECK(!hasIndex(R,ci));
        CHECK(hasIndex(R,L1));
        CHECK(hasIndex(R,L2));
        CHECK(hasIndex(R,S1));
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.dim(); ++i1)
        for(int i2 = 1; i2 <= L2.dim(); ++i2)
        for(int j1 = 1; j1 <= S1.dim(); ++j1)
            {
            CHECK_CLOSE( T.elt(L1(i1),L2(i2),S1(j1)), R.elt(L1(i1),L2(i2),S1(j1)) );
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
        for(int i1 = 1; i1 <= L1.dim(); ++i1)
        for(int i2 = 1; i2 <= L2.dim(); ++i2)
        for(int j1 = 1; j1 <= S1.dim(); ++j1)
        for(int j2 = 1; j2 <= S2.dim(); ++j2)
            {
            CHECK_CLOSE( T.elt(L1(i1),L2(i2),S1(j1),S2(j2)), R.elt(L1(i1),L2(i2),S1(j1),S2(j2)) );
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
        for(int i1 = 1; i1 <= L1.dim(); ++i1)
        for(int i2 = 1; i2 <= L2.dim(); ++i2)
        for(int j1 = 1; j1 <= S1.dim(); ++j1)
        for(int j2 = 1; j2 <= S2.dim(); ++j2)
            {
            CHECK_CLOSE( T.elt(L1(i1),L2(i2),S1(j1),S2(j2)), R.elt(L1(i1),L2(i2),S1(j1),S2(j2)) );
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
        CHECK(!hasIndex(R,L1));
        CHECK(!hasIndex(R,S1));
        CHECK(!hasIndex(R,L2));
        CHECK(hasIndex(R,S2));
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R = dag(C)*R; //uncombine
        CHECK(!hasIndex(R,ci));
        CHECK(hasIndex(R,L1));
        CHECK(hasIndex(R,L2));
        CHECK(hasIndex(R,S1));
        CHECK(hasIndex(R,S2));
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.dim(); ++i1)
        for(int i2 = 1; i2 <= L2.dim(); ++i2)
        for(int j1 = 1; j1 <= S1.dim(); ++j1)
        for(int j2 = 1; j2 <= S2.dim(); ++j2)
            {
            CHECK_CLOSE( T.elt(L1(i1),L2(i2),S1(j1),S2(j2)), R.elt(L1(i1),L2(i2),S1(j1),S2(j2)) );
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
        CHECK(hasIndex(R,L1));
        CHECK(hasIndex(R,L2));
        CHECK(!hasIndex(R,S1));
        CHECK(!hasIndex(R,S2));
        CHECK_CLOSE(norm(T),norm(R));
        CHECK(div(T) == div(R));

        R *= dag(C); //uncombine
        //Check that R equals original T
        for(int i1 = 1; i1 <= L1.dim(); ++i1)
        for(int i2 = 1; i2 <= L2.dim(); ++i2)
        for(int j1 = 1; j1 <= S1.dim(); ++j1)
        for(int j2 = 1; j2 <= S2.dim(); ++j2)
            {
            CHECK_CLOSE( T.elt(L1(i1),L2(i2),S1(j1),S2(j2)), R.elt(L1(i1),L2(i2),S1(j1),S2(j2)) );
            }
        }

    SECTION("Combine / Uncombined 7")
        {
        auto T = randomTensor(QN(),L1,S1,L2,S2);
        auto C = combiner(L2,L1);
        auto R = T*C;
        CHECK_CLOSE(norm(R),norm(T));

        auto U = dag(C)*R;
        for(auto l1 : range1(L1.dim()))
        for(auto l2 : range1(L1.dim()))
        for(auto s1 : range1(S1.dim()))
        for(auto s2 : range1(S2.dim()))
            {
            CHECK_CLOSE(U.elt(L1(l1),L2(l2),S1(s1),S2(s2)),T.elt(L1(l1),L2(l2),S1(s1),S2(s2)));
            }
        }

    SECTION("Fragmented IQIndex Combiner Test 1")
        {
        auto i1 = IQIndex(Index(2),QN(+1),
                          Index(2),QN(+0),
                          Index(2),QN(-1),
                          Index(2),QN(+2),
                          Index(2),QN(-2),
                          Index(2),QN(+0),
                          Index(2),QN(-1));

        auto i2 = IQIndex(Index(2),QN(+1),
                          Index(4),QN(-1),
                          Index(2),QN(+0),
                          Index(3),QN(+2),
                          Index(2),QN(+0),
                          Index(2),QN(-1),
                          Index(3),QN(+1));

        auto i3 = IQIndex(Index(3),QN(-1),
                          Index(2),QN(+0),
                          Index(2),QN(+1),
                          Index(2),QN(+0),
                          Index(4),QN(-1),
                          Index(2),QN(+0),
                          Index(2),QN(+1));


        SECTION("Test 1")
            {
            auto flux = QN(-2);
            auto T = randomTensor(flux,i1,prime(i3),i2,i3,prime(i2));
            auto C = combiner(i1,i2,prime(i2));
            auto R = C * T;
            auto ci = commonIndex(R,C);
            
            ////Uncombine
            auto nT = dag(C) * R;
            CHECK(div(T) == div(R));
            CHECK(div(T) == div(nT));
            CHECK(norm(T-nT) < 1E-11);
            }
        SECTION("Test 2")
            {
            auto flux = QN(-2);
            auto T = randomTensor(flux,i1,prime(i3),i2,i3,prime(i2));
            auto C = combiner(i3,i1);
            auto R = C * T;

            ////Uncombine
            auto nT = dag(C) * R;

            CHECK(div(T) == div(R));
            CHECK(div(T) == div(nT));
            CHECK(norm(T-nT) < 1E-11);
            }

        SECTION("Combine Twice and Uncombine")
            {
            auto T = randomTensor(QN(),i1,i2,prime(i2),prime(i1),i3);
            auto C1 = combiner(i1,i3);
            auto C2 = combiner(i2,prime(i1));

            auto cT = C1 * T * C2;
            auto ci1 = commonIndex(C1,cT);
            auto ci2 = commonIndex(C2,cT);
            CHECK(ci1.dim() == (i1.dim()*i3.dim()));
            CHECK(ci2.dim() == (i2.dim()*i1.dim()));
            CHECK(hasIndex(cT,prime(i2)));
            CHECK(div(T) == div(cT));
            CHECK(std::fabs(norm(cT)-norm(T)) < 1E-11);

            auto uT = dag(C1)*cT*dag(C2);
            CHECK(div(T) == div(uT));
            CHECK(norm(uT-T) < 1E-11);
            }
        }

    SECTION("Combiner Arrow Regression Test")
        {
        IQIndex s9(Index{1},QN(+1),
                   Index{1},QN(-1));

        IQIndex hl8(Index{3},QN(0),
                    Index{1},QN(-2),
                    Index{1},QN(+2));

        IQIndex L(Index{1},QN(0));

        auto A = randomTensor(QN(),dag(s9),prime(s9),dag(hl8),dag(L));
        auto C = combiner(dag(s9),prime(s9),dag(L));

        auto R = C * A;

        CHECK(fabs(norm(R)-norm(A)) < 1E-11);

        auto nA1 = dag(C) * R;
        CHECK(norm(nA1-A) < 1E-11);

        auto nA2 = R * dag(C);
        CHECK(norm(nA2-A) < 1E-11);

        }

    SECTION("Combiner Single Index Regression Test")
        {
        auto s1 = IQIndex(Index(1),QN("Sz=", 0,"Pf=",0),
                          Index(1),QN("Sz=",+1,"Pf=",1),
                          Index(1),QN("Sz=",-1,"Pf=",1),
                          Index(1),QN("Sz=", 0,"Pf=",0));
        auto s2 = IQIndex(Index(1),QN("Sz=", 0,"Pf=",0),
                          Index(1),QN("Sz=",+1,"Pf=",1),
                          Index(1),QN("Sz=",-1,"Pf=",1),
                          Index(1),QN("Sz=", 0,"Pf=",0));
        auto T = randomTensor(QN{},s1,s2);
        auto cmb = combiner(s1);
        auto ci = cmb.inds().front();
        auto Tc = cmb*T;
        //Following line was causing an error:
        Tc *= dag(prime(Tc,ci));
        }


    } //Combiner

SECTION("Diag IQTensor Contraction")
{
SECTION("Diag All Same")
    {
    auto d = delta(dag(S1),S2); //all diag elements same
    CHECK(typeOf(d) == QType::QDiagRealAllSame);

    auto T = randomTensor(QN(),S1,prime(S1,2));
    auto R = d*T;
    CHECK(hasIndex(R,S2));
    CHECK(hasIndex(R,prime(S1,2)));

    for(auto j1 : range1(S1))
    for(auto j2 : range1(S2))
        {
        CHECK_CLOSE(R.elt(prime(S1,2)(j1),S2(j2)), T.elt(prime(S1,2)(j1),S1(j2)));
        }
    }

SECTION("Trace")
    {
    SECTION("Case 1")
        {
        auto T = randomTensor(QN(),S1,dag(S2),S3,S4);
        auto d = delta(dag(S1),S2);
        auto R = d*T;

        for(auto i3 : range1(S3))
        for(auto i4 : range1(S4))
            {
            Real val = 0;
            for(auto i12 : range1(S1))
                {
                val += T.elt(S1(i12),S2(i12),S3(i3),S4(i4));
                }
            CHECK_CLOSE(val,R.elt(S3(i3),S4(i4)));
            }
        }

    SECTION("Case 2")
        {
        auto ti = IQIndex(Index(2),QN(0));

        auto T = randomTensor(QN(),S1,dag(S3),dag(ti));
        auto d = delta(S3,ti,dag(S1));
        auto R = d*T;

        Real val = 0;
        for(auto j : range1(ti))
            {
            val += T.elt(S1(j),S3(j),ti(j));
            }
        CHECK_CLOSE(val,R.elt());
        }
    }

SECTION("Tie Indices with Diag IQTensor")
    {
    SECTION("Case 1")
        {
        //both tied indices have In arrows
        auto T = randomTensor(QN(),S1,S2,S3,S4);
        auto ti = IQIndex(Index(1),QN(-2),
                          Index(1),QN(+2));
        auto tie = delta(dag(S1),dag(S2),ti);
        auto R = T*tie;

        for(auto i : range1(ti))
        for(auto j3 : range1(S3))
        for(auto j4 : range1(S4))
            {
            CHECK_CLOSE(T.elt(S1(i),S2(i),S3(j3),S4(j4)), R.elt(ti(i),S3(j3),S4(j4)));
            }
        }

    SECTION("Case 2")
        {
        //tied indices have mixed arrows
        auto T = randomTensor(QN(),S1,S2,dag(S3),S4);
        auto ti = IQIndex(Index(2),QN(0));
        auto tie = delta(dag(S1),S3,ti);

        auto R = T*tie;

        for(auto i : range1(ti))
        for(auto j2 : range1(S2))
        for(auto j4 : range1(S4))
            {
            CHECK_CLOSE(T.elt(S1(i),S2(j2),S3(i),S4(j4)), R.elt(ti(i),S2(j2),S4(j4)));
            }
        }

    SECTION("Case 3")
        {
        //tied indices have mixed arrows
        //similar to above case but different index
        //order for tie IQTensor
        auto T = randomTensor(QN(),S1,S2,dag(S3),S4);
        auto ti = IQIndex(Index(2),QN(0));
        auto tie = delta(ti,dag(S1),S3);

        auto R = T*tie;

        for(auto i : range1(ti))
        for(auto j2 : range1(S2))
        for(auto j4 : range1(S4))
            {
            CHECK_CLOSE(T.elt(S1(i),S2(j2),S3(i),S4(j4)), R.elt(ti(i),S2(j2),S4(j4)));
            }
        }

    SECTION("Case 4")
        {
        //tied indices have mixed arrows
        //similar to above case but yet another
        //index order for tie IQTensor
        auto T = randomTensor(QN(),S1,S2,dag(S3),S4);
        auto ti = IQIndex(Index(2),QN(0));
        auto tie = delta(S3,ti,dag(S1));

        auto R = T*tie;

        for(auto i : range1(ti))
        for(auto j2 : range1(S2))
        for(auto j4 : range1(S4))
            {
            CHECK_CLOSE(T.elt(S1(i),S2(j2),S3(i),S4(j4)), R.elt(ti(i),S2(j2),S4(j4)));
            }
        }
    }


//SECTION("Contract All Dense Inds; Diag result")
//    {
//    auto T = randomTensor(QN{},L1,prime(L1));
//    
//    auto d = delta(dag(L1),dag(prime(L1)),prime(L1,2),prime(L1,3));
//    auto R = d*T;
//    //CHECK(typeOf(R) == Type::DiagReal);
//    //CHECK(hasIndex(R,L));
//    //auto minjkl = std::min(std::min(J.dim(),K.dim()),L.dim());
//    //for(long j = 1; j <= minjkl; ++j)
//    //    CHECK_CLOSE(R.elt(L(j)), T.elt(J(j),K(j)));
//    }

//SECTION("Contract All Dense Inds; Diag Scalar result")
//    {
//    auto T = randomTensor(J,K);
//
//    auto d1 = delta(J,K);
//    auto R = d1*T;
//    CHECK(typeOf(R) == Type::DiagRealAllSame);
//    Real val = 0;
//    auto minjk = std::min(J.dim(),K.dim());
//    for(long j = 1; j <= minjk; ++j)
//        val += T.elt(J(j),K(j));
//    CHECK_CLOSE(R.elt(),val);
//
//    auto data = randomData(minjk);
//    auto d2 = diagTensor(data,J,K);
//    R = d2*T;
//    CHECK(typeOf(R) == Type::DiagRealAllSame);
//    val = 0;
//    for(long j = 1; j <= minjk; ++j)
//        val += data.at(j-1)*T.elt(J(j),K(j));
//    CHECK_CLOSE(R.elt(),val);
//    }

SECTION("Single IQIndex Replacement")
    {
    auto d = delta(dag(S1),S2);
    CHECK(isIQTensor(d));
    auto T = randomTensor(QN{},S1,prime(S1));
    CHECK(isIQTensor(T));

    auto R = d*T;
    CHECK(isIQTensor(R));
    CHECK(hasIndex(R,S2));
    CHECK(hasIndex(R,prime(S1)));

    for(auto i1 : range1(S1.dim()))
    for(auto i2 : range1(S1.dim()))
        {
        CHECK(T.elt(S1(i1),prime(S1)(i2)) == R.elt(S2(i1),prime(S1)(i2)));
        }
    }

SECTION("Replace two IQIndex's with three")
    {
    auto ti = IQIndex(Index(2),QN(0));
    auto T = randomTensor(QN(),S1,dag(S3));
    auto d = delta(S3,dag(S1),S4,dag(S2),ti);
    auto R = d*T;
    CHECK(typeOf(R) == QType::QDiagReal);

    for(auto j : range1(ti))
        {
        CHECK_CLOSE(R.elt(S2(j),S4(j),ti(j)), T.elt(S1(j),S3(j)));
        }
    }
}


SECTION("Scalar")
    {
    auto T1 = randomTensor(QN(),L1,L2,S1,S2);
    auto T2 = randomTensor(QN(),S1,L2,S2,L1);
    auto S = T1*dag(T2);
    CHECK(S.r() == 0);

    Real val = 0;
    for(int i1 = 1; i1 <= L1.dim(); ++i1)
    for(int i2 = 1; i2 <= L2.dim(); ++i2)
    for(int j1 = 1; j1 <= S1.dim(); ++j1)
    for(int j2 = 1; j2 <= S2.dim(); ++j2)
        {
        val += T1.elt(L1(i1),L2(i2),S1(j1),S2(j2))*T2.elt(L1(i1),L2(i2),S1(j1),S2(j2));
        }

    CHECK_CLOSE(val, S.elt());
    CHECK_CLOSE(fabs(val), norm(S));
    }

SECTION("Contraction with Scalar")
    {
    auto T1 = randomTensor(QN(),L1,L2,S1,S2);
    auto T2 = randomTensor(QN(),S1,L2,S2,L1);
    auto S = T1*dag(T2);
    CHECK(S.r() == 0);
    auto T3 = S * T1;

    auto z = S.eltC();
    for(int i1 = 1; i1 <= L1.dim(); ++i1)
    for(int i2 = 1; i2 <= L2.dim(); ++i2)
    for(int j1 = 1; j1 <= S1.dim(); ++j1)
    for(int j2 = 1; j2 <= S2.dim(); ++j2)
        {
        auto val = z*T1.elt(L1(i1),L2(i2),S1(j1),S2(j2));
        CHECK_CLOSE(val,T3.elt(L1(i1),L2(i2),S1(j1),S2(j2)));
        }

    //Make real scalar tensor from constructor
    auto Q = IQTensor(2);
    auto q = Q.elt();
    auto T4 = Q * T1;
    for(int i1 = 1; i1 <= L1.dim(); ++i1)
    for(int i2 = 1; i2 <= L2.dim(); ++i2)
    for(int j1 = 1; j1 <= S1.dim(); ++j1)
    for(int j2 = 1; j2 <= S2.dim(); ++j2)
        {
        auto val = q*T1.elt(L1(i1),L2(i2),S1(j1),S2(j2));
        CHECK_CLOSE(val,T4.elt(L1(i1),L2(i2),S1(j1),S2(j2)));
        }

    ////Make complex scalar tensor from constructor
    //auto Q = IQTensor(1.+3_i);
    //auto q = Q.eltC();
    //auto T4 = Q * T1;
    //for(int i1 = 1; i1 <= L1.dim(); ++i1)
    //for(int i2 = 1; i2 <= L2.dim(); ++i2)
    //for(int j1 = 1; j1 <= S1.dim(); ++j1)
    //for(int j2 = 1; j2 <= S2.dim(); ++j2)
    //    {
    //    auto val = q*T1.eltC(L1(i1),L2(i2),S1(j1),S2(j2));
    //    CHECK_CLOSE(val,T4.eltC(L1(i1),L2(i2),S1(j1),S2(j2)));
    //    }
    }

SECTION("Scalar Storage")
    {
    auto S1 = IQTensor(1.);
    CHECK_CLOSE(S1.elt(),1.);

    auto S2 = IQTensor(1.)*2.;
    CHECK_CLOSE(S2.elt(),2.);

    auto ZA = IQTensor(1._i);
    CHECK_CLOSE(ZA.eltC(),1._i);

    auto ZB = IQTensor(-1.+2._i);
    CHECK_CLOSE(ZB.eltC(),-1+2._i);

    SECTION("Set")
        {
        S1.set(4.5);
        CHECK_CLOSE(S1.elt(),4.5);
        S1.set(1.+3._i);
        CHECK(isComplex(S1));
        CHECK_CLOSE(S1.eltC(),1.+3._i);

        ZA.set(2.-3._i);
        CHECK_CLOSE(ZA.eltC(),2.-3._i);
        ZA.set(3.0);
        CHECK(isReal(ZA));
        CHECK_CLOSE(ZA.elt(),3.);
        }

    SECTION("Norm")
        {
        CHECK_CLOSE(norm(S1),1.);
        CHECK_CLOSE(norm(S2),2.);
        auto Sn2 = IQTensor(-2.);
        CHECK_CLOSE(norm(Sn2),2.);

        CHECK_CLOSE(norm(ZA),std::norm(ZA.eltC()));
        }

    auto T = randomTensor(QN(),L1,L2);
    auto TC = randomTensorC(QN(),L1,L2);

    SECTION("Multiply on right")
        {
        auto R = T*S1;
        CHECK(norm(R-T) < 1E-12);
        R = T*S2;
        CHECK(norm(R-2*T) < 1E-12);
        R = TC*S2;
        CHECK(norm(R-2*TC) < 1E-12);

        R = T*ZA;
        CHECK(isComplex(R));
        CHECK(norm(R-ZA.eltC()*T) < 1E-12);
        R = TC*ZA;
        CHECK(norm(R-ZA.eltC()*TC) < 1E-12);
        }
    SECTION("Multiply on left")
        {
        auto R = S1*T;
        CHECK(norm(R-T) < 1E-12);
        R = S2*T;
        CHECK(norm(R-2*T) < 1E-12);
        R = S2*TC;
        CHECK(norm(R-2*TC) < 1E-12);

        R = ZA*T;
        CHECK(isComplex(R));
        CHECK(norm(R-ZA.eltC()*T) < 1E-12);

        R = ZA*TC;
        CHECK(norm(R-ZA.eltC()*TC) < 1E-12);
        }
    SECTION("Add & Subtract")
        {
        auto R = S1 + S2;
        CHECK_CLOSE(R.elt(),3.);

        R = ZA+ZB;
        CHECK_CLOSE(R.eltC(),ZA.eltC()+ZB.eltC());

        R = S1 - S2;
        CHECK_CLOSE(R.elt(),-1.);

        R = ZA-ZB;
        CHECK_CLOSE(R.eltC(),ZA.eltC()-ZB.eltC());
        }
    }

SECTION("Mixed Storage")
    {
    auto s = IQIndex(Index(1),QN(+1),
                     Index(1),QN(-1));
    auto T = mixedIQTensor(s,prime(s));
    T.set(s(1),prime(s)(1),11);
    T.set(s(1),prime(s)(2),12);
    T.set(s(2),prime(s)(1),21);
    T.set(s(2),prime(s)(2),22);

    auto t = ITensor(T);
    for(auto i : range1(s))
    for(auto j : range1(s))
        {
        CHECK_CLOSE(t.elt(s(i),prime(s)(j)),T.elt(s(i),prime(s)(j)));
        }
    }

SECTION("isEmpty Function")
    {
    IQTensor T;
    CHECK(isEmpty(T));

    T = randomTensor(QN{},S1,S2,S3,S4);
    CHECK(not isEmpty(T));

    auto D = delta(S3,dag(S1),S4,dag(S2));
    CHECK(typeOf(D) == QType::QDiagRealAllSame);
    CHECK(not isEmpty(D));
    }

SECTION("Permute Test")
{
auto i = IQIndex(Index(1),QN(-1),
                 Index(1),QN(+1));
auto j = IQIndex(Index(2),QN(-1),
                 Index(2),QN( 0),
                 Index(2),QN(+1));
auto k = IQIndex(Index(2),QN(-1),
                 Index(3),QN( 0),
                 Index(2),QN(+1));
auto jp = prime(j);

auto IT = randomTensor(QN(0),i,j,jp,dag(k));

auto O1 = permute(IT,jp,k,j,i);
CHECK(IT.inds().index(1)==O1.inds().index(4));
CHECK(IT.inds().index(2)==O1.inds().index(3));
CHECK(IT.inds().index(3)==O1.inds().index(1));
CHECK(IT.inds().index(4)==O1.inds().index(2));
CHECK(IT.inds().index(1).dir()==O1.inds().index(4).dir());
CHECK(IT.inds().index(2).dir()==O1.inds().index(3).dir());
CHECK(IT.inds().index(3).dir()==O1.inds().index(1).dir());
CHECK(IT.inds().index(4).dir()==O1.inds().index(2).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(IT.elt(i(ii),j(jj),jp(jjp),k(kk)),O1.elt(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O2 = permute(IT,j,i,k,jp);
CHECK(IT.inds().index(1)==O2.inds().index(2));
CHECK(IT.inds().index(2)==O2.inds().index(1));
CHECK(IT.inds().index(3)==O2.inds().index(4));
CHECK(IT.inds().index(4)==O2.inds().index(3));
CHECK(IT.inds().index(1).dir()==O2.inds().index(2).dir());
CHECK(IT.inds().index(2).dir()==O2.inds().index(1).dir());
CHECK(IT.inds().index(3).dir()==O2.inds().index(4).dir());
CHECK(IT.inds().index(4).dir()==O2.inds().index(3).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(IT.elt(i(ii),j(jj),jp(jjp),k(kk)),O2.elt(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto CIT = randomTensorC(QN(0),i,j,jp,k);

auto O3 = permute(CIT,jp,k,i,j);
CHECK(CIT.inds().index(1)==O3.inds().index(3));
CHECK(CIT.inds().index(2)==O3.inds().index(4));
CHECK(CIT.inds().index(3)==O3.inds().index(1));
CHECK(CIT.inds().index(4)==O3.inds().index(2));
CHECK(CIT.inds().index(1).dir()==O3.inds().index(3).dir());
CHECK(CIT.inds().index(2).dir()==O3.inds().index(4).dir());
CHECK(CIT.inds().index(3).dir()==O3.inds().index(1).dir());
CHECK(CIT.inds().index(4).dir()==O3.inds().index(2).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(CIT.eltC(i(ii),j(jj),jp(jjp),k(kk)),O3.eltC(i(ii),j(jj),jp(jjp),k(kk)));
    }

}

SECTION("Permute Test: Dots Syntax")
{
auto i = IQIndex(Index(1),QN(-1),
                 Index(1),QN(+1));
auto j = IQIndex(Index(2),QN(-1),
                 Index(2),QN( 0),
                 Index(2),QN(+1));
auto k = IQIndex(Index(2),QN(-1),
                 Index(3),QN( 0),
                 Index(2),QN(+1));
auto jp = prime(j);

auto IT = randomTensor(QN(0),i,j,jp,dag(k));

auto O1 = permute(IT,"...",i);
CHECK(IT.inds().index(1)==O1.inds().index(4));
CHECK(IT.inds().index(1).dir()==O1.inds().index(4).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(IT.elt(i(ii),j(jj),jp(jjp),k(kk)),O1.elt(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O2 = permute(IT,"...",j,i);
CHECK(IT.inds().index(1)==O2.inds().index(4));
CHECK(IT.inds().index(2)==O2.inds().index(3));
CHECK(IT.inds().index(1).dir()==O2.inds().index(4).dir());
CHECK(IT.inds().index(2).dir()==O2.inds().index(3).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(IT.elt(i(ii),j(jj),jp(jjp),k(kk)),O2.elt(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O3 = permute(IT,"...",k,j,i);
CHECK(IT.inds().index(1)==O3.inds().index(4));
CHECK(IT.inds().index(2)==O3.inds().index(3));
CHECK(IT.inds().index(3)==O3.inds().index(1));
CHECK(IT.inds().index(4)==O3.inds().index(2));
CHECK(IT.inds().index(1).dir()==O3.inds().index(4).dir());
CHECK(IT.inds().index(2).dir()==O3.inds().index(3).dir());
CHECK(IT.inds().index(3).dir()==O3.inds().index(1).dir());
CHECK(IT.inds().index(4).dir()==O3.inds().index(2).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(IT.elt(i(ii),j(jj),jp(jjp),k(kk)),O3.elt(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O4 = permute(IT,k,"...");
CHECK(IT.inds().index(4)==O4.inds().index(1));
CHECK(IT.inds().index(4).dir()==O4.inds().index(1).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(IT.elt(i(ii),j(jj),jp(jjp),k(kk)),O4.elt(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O5 = permute(IT,k,jp,i,"...");
CHECK(IT.inds().index(1)==O5.inds().index(3));
CHECK(IT.inds().index(2)==O5.inds().index(4));
CHECK(IT.inds().index(3)==O5.inds().index(2));
CHECK(IT.inds().index(4)==O5.inds().index(1));
CHECK(IT.inds().index(1).dir()==O5.inds().index(3).dir());
CHECK(IT.inds().index(2).dir()==O5.inds().index(4).dir());
CHECK(IT.inds().index(3).dir()==O5.inds().index(2).dir());
CHECK(IT.inds().index(4).dir()==O5.inds().index(1).dir());
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(IT.elt(i(ii),j(jj),jp(jjp),k(kk)),O5.elt(i(ii),j(jj),jp(jjp),k(kk)));
    }

}

SECTION("Set and Get with int")
{
auto I = IQIndex(Index(1),QN(+1),
                 Index(1),QN(-1));
auto J = IQIndex(Index(2),QN(+1),
                 Index(2),QN(-1));
auto T = IQTensor(I,J);
T.set(2,1,21);
CHECK_CLOSE(T.elt(J(1),I(2)),21);
CHECK_CLOSE(T.elt(2,1),21);
}

SECTION("Set and Get with long int")
{
auto I = IQIndex(Index(1),QN(+1),
                 Index(1),QN(-1));
auto J = IQIndex(Index(2),QN(+1),
                 Index(2),QN(-1));
auto T = IQTensor(I,J);
long int i1 = 1,
         i2 = 2;
T.set(i2,i1,21);
CHECK_CLOSE(T.elt(J(1),I(2)),21);
CHECK_CLOSE(T.elt(i2,i1),21);
}

/*
//TODO: this function should not ignore primes and tags
SECTION("Reindex")
    {

    auto T = randomTensor(QN(),dag(S1),prime(S1),dag(S2),prime(S2));

    auto nT = reindex(T,S1,S3,S2,S4);
    
    CHECK(order(nT) == 4);
    CHECK(hasIndex(nT,S3));
    CHECK(hasIndex(nT,prime(S3)));
    CHECK(hasIndex(nT,S4));
    CHECK(hasIndex(nT,prime(S4)));

    auto J4 = IQIndex(Index(1,"J4,Site"),QN(-1),
                      Index(1,"J4,Site"),QN(0),
                      Index(1,"J4,Site"),QN(+1));

    CHECK_THROWS(nT = reindex(T,S1,S3,S2,J4));
    }
*/

//SECTION("Non-contracting product")
//    {
//    SECTION("Case 1")
//        {
//        //This use case of ncprod may not
//        //lead to a C with a well defined divergence!
//        
//        auto s = IQIndex("s",Index{"Up s",2,Site},QN(+1),
//                             Index{"Dn s",2,Site},QN(-1));
//        auto t = IQIndex("t",Index{"Up t",2,Site},QN(+1),
//                             Index{"Dn t",2,Site},QN(-1));
//        auto u = IQIndex("u",Index{"Up u",2,Site},QN(+1),
//                             Index{"Dn u",2,Site},QN(-1));
//
//        auto h = IQIndex("h",Index{"h 0",3},QN(0),
//                             Index{"h-1",3},QN(-1),
//                             Index{"h-2",3},QN(-2),
//                             Index{"h+2",3},QN(+2));
//
//        auto k = IQIndex("k",Index{"k+1",3},QN(1),
//                             Index{"k-1",3},QN(-1),
//                             Index{"k+2",3},QN(+2));
//
//        auto A = randomTensor(QN{-1},s,h);
//        auto B = randomTensor(QN{0},h,k,t);
//
//        auto C = A/B;
//
//        auto diff = 0.;
//        for(auto S : range1(s.dim()))
//        for(auto H : range1(h.dim()))
//        for(auto K : range1(k.dim()))
//        for(auto T : range1(t.dim()))
//            {
//            diff += C.elt(t(T),s(S),h(H),k(K)) - A.elt(s(S),h(H))*B.elt(h(H),k(K),t(T));
//            }
//        CHECK(diff < 1E-13);
//        }
//
//    }



//SECTION("TieIndices")
//    {
//    IQTensor D1 = tieIndices(D,L1,prime(L1),L1);
//
//    for(int k1 = 1; k1 <= L1.dim(); ++k1)
//    for(int k2 = 1; k2 <= L1.dim(); ++k2)
//    for(int k3 = 1; k3 <= S1.dim(); ++k3)
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
//    for(int j2 = 1; j2 <= S1.dim(); ++j2)
//    for(int j1 = 1; j1 <= L1.dim(); ++j1)
//        {
//        Real val = 0;
//        for(int k1 = 1; k1 <= L1.dim(); ++k1)
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
