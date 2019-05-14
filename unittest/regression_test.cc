#include "test.h"
#include "itensor/decomp.h"

using namespace itensor;
using namespace std;

//
//The tests in this suite are
//regression tests, meaning they
//are designed to check for a specific
//bug that came up in the use of the library.
//Keeping these tests around is supposed
//to prevent the same bugs from coming up again.
//

//TODO: doTask not defined for task Contract and storage types QDenseReal DenseReal
TEST_CASE("ITensor Times IndexVal (QN)")
    {
    Index s(QN(+1),1,
            QN(-1),1);

    Index l(4);
    ITensor T(l);
    T.randomize();

    ITensor R = T * setElt(s(2));

    REQUIRE(hasIndex(R,s));
    CHECK(R.elt(l(1),s(1)) == 0);
    CHECK(R.elt(l(2),s(1)) == 0);
    CHECK(R.elt(l(3),s(1)) == 0);
    CHECK(R.elt(l(4),s(1)) == 0);
    }

TEST_CASE("ITensor from IndexVal (QN)")
    {
    Index s(QN(+1),1,
            QN(-1),1);

    auto T1 = setElt(s(1));
    CHECK(T1.elt(s(1)) == 1);
    CHECK(T1.elt(s(2)) == 0);

    auto T2 = setElt(s(2));
    CHECK(T2.elt(s(1)) == 0);
    CHECK(T2.elt(s(2)) == 1);
    }

TEST_CASE("CombinerOrder")
    {
    Index a(2),c(2);

    ITensor U(a,c);
    U.set(a(1),c(2),1);
    //PrintDat(U);

    auto [C,ci] = combiner(c,a);
    auto CU = C * U;
    CHECK(hasIndex(CU,ci));
    auto UU = C * CU;

    CHECK(norm(U-UU) < 1E-10);

    auto [D,di] = combiner(a,c);
    CU = D * U;
    CHECK(hasIndex(CU,di));
    UU = D * CU;

    CHECK(norm(U-UU) < 1E-10);
    }

TEST_CASE("SVDIndexOrder")
    {
    Index a(2),b(1),c(2);

    ITensor z(c,a,b);
    z.set(c(1),a(2),b(1),1);
    z.set(c(1),a(1),b(1),1);


    ITensor V(b);

    //ITensor U(c,a); //This order succeeds
    ITensor U(a,c); //This order was failing

    ITensor D;
    //Globals::debug2() = true;
    svd(z,U,D,V);
    //Globals::debug2() = false;

    //PrintDat(U);

    //Combiner c1(a,c),c2(c,a);

    //U = c2*U;
    //PrintDat(U);
    //U = c2*U;
    //PrintDat(U);

    ITensor nz = U*D*V;
    CHECK(norm(z-nz) < 1E-10);

    }

//TODO: this doesn't work anymore now that IQIndex is not
//      made up of Indices
//TEST_CASE("SVDArrows")
//    {
//    Index L(QN(1,1),2,In),R(QN(1,1),2,Out);
//
//    ITensor AA(L,R);
//
//    ITensor block(l,r);
//    randomize(block);
//    AA += block;
//
//    const QN Zero;
//    CHECK_EQUAL(div(AA),Zero);
//
//    ITensor U(L),V(R);
//    ITensor D;
//    svd(AA,U,D,V);
//
//    CHECK_EQUAL(div(U),Zero);
//    CHECK_EQUAL(div(V),Zero);
//    }

//TEST_CASE("ExpandIndex")
//    {
//    //
//    //ITensor::expandIndex was
//    //failing when the first Index 
//    //had m == 1 and the second or
//    //third Index was to be expanded.
//    //The reason is that expandIndex
//    //assumed a stable Index order
//    //whereas the ITensor constructor
//    //moves m==1's to the back
//    //
//    Index l("l");
//
//    Index emp("emp"),occ("occ");
//    Index S("S",emp,QN(0,0),
//                occ,QN(1,0),
//                Out);
//
//    ITensor oo(l,occ,prime(occ));
//    oo(l(1),occ(1),prime(occ)(1)) = 1;
//
//    oo.expandIndex(occ,S,offset(S,occ));
//
//    oo.expandIndex(prime(occ),prime(S),offset(S,occ));
//
//    CHECK_CLOSE(0,oo(S(1),prime(S)(1),l(1)),1E-5);
//    CHECK_CLOSE(1,oo(S(2),prime(S)(2),l(1)),1E-5);
//    }

//TODO: create toDense(ITensor) function
TEST_CASE("ConvertQDenseToDenseITensor")
    {
    Index L(QN(),1,Out);
    Index emp(1),occ(1);
    Index S(QN(0),1,
            QN(1),1,
            Out);

    auto T = ITensor(L,dag(S),prime(S));
    T.set(L(1),S(1),prime(S)(1),1);
    T.set(L(1),S(2),prime(S)(2),1);

    auto t = removeQNs(T);

    CHECK_CLOSE(1,t.elt(L(1),S(1),prime(S)(1)));
    CHECK_CLOSE(1,t.elt(L(1),S(2),prime(S)(2)));
    CHECK_CLOSE(0,t.elt(L(1),S(2),prime(S)(1)));
    CHECK_CLOSE(0,t.elt(L(1),S(2),prime(S)(1)));
    }

//TEST_CASE("IndexOrder")
//    {
//    //
//    //The ITensor contracting product code,
//    //specifically the toMatrixProd part,
//    //was mis-identifying certain products 
//    //as matrix products when they were not.
//    //
//
//    //Globals::debug1() = true;
//
//    Index l("l",2),r("r",2),s("s",2,Site);
//
//    ITensor HP(l,s,r);
//    HP(l(1),s(2),r(1)) = 4.30425;
//
//    ITensor phi(l,s,r);
//    phi(l(1),s(2),r(1)) = -0.341723;
//
//    ITensor phialt(l,r,s);
//    phialt.assignFrom(phi);
//
//    //PrintDat(phi);
//    //PrintDat(phialt);
//
//    CHECK_CLOSE((phi-phialt).norm(),0,1E-3);
//
//    ITensor order2 = HP * phi;
//    //Print(order2.val0());
//
//    //cout << endl << endl;
//
//    ITensor order2alt = HP * phialt;
//    //Print(order2alt.val0());
//
//
//    CHECK_CLOSE(order2.val0(),order2alt.val0(),1E-5);
//
//    }

TEST_CASE("ComplexAddition")
    {
    //EMS Oct 21, 2013
    //Bug was happening because a below has different
    //Index order from b but complex addition code 
    //did not account for this!
    Index L1(6),
          S1(2);

    ITensor a(L1,S1);
    a.set(L1(1),S1(2),21);
    a.set(L1(2),S1(2),22);
    a.scaleTo(-1.2243);

    ITensor ca = a*Complex_i;

    ITensor b(S1,L1);
    b.set(S1(1),L1(3),13);
    b.set(S1(1),L1(4),14);
    b.scaleTo(9.3435);

    ITensor r = b + ca;

    CHECK(norm(realPart(r)-b) < 1E-12);
    CHECK(norm(imagPart(r)-a) < 1E-12);
    }

