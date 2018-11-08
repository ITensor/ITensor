#include "test.h"
#include "itensor/decomp.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using namespace std;

TEST_CASE("Additional SVD implementation Tests")
{

SECTION("dgesdd test")
    {
    Index K  = Index("k",5);
    Index Kp = prime(K);

    ITensor A(K,Kp);
    A.set(K(1),Kp(1), 6.80);
    A.set(K(1),Kp(2),-6.05);
    A.set(K(1),Kp(3),-0.45);
    A.set(K(1),Kp(4), 8.32);
    A.set(K(1),Kp(5),-9.67); 

    A.set(K(2),Kp(1),-2.11);
    A.set(K(2),Kp(2),-3.30);
    A.set(K(2),Kp(3), 2.58);
    A.set(K(2),Kp(4), 2.71);
    A.set(K(2),Kp(5),-5.14); 

    A.set(K(3),Kp(1), 5.66);
    A.set(K(3),Kp(2), 5.36);
    A.set(K(3),Kp(3),-2.70);
    A.set(K(3),Kp(4), 4.35);
    A.set(K(3),Kp(5),-7.26); 

    A.set(K(4),Kp(1), 5.97);
    A.set(K(4),Kp(2),-4.44);
    A.set(K(4),Kp(3), 0.27);
    A.set(K(4),Kp(4),-7.17);
    A.set(K(4),Kp(5), 6.08); 

    A.set(K(5),Kp(1), 8.23);
    A.set(K(5),Kp(2), 1.08);
    A.set(K(5),Kp(3), 9.04);
    A.set(K(5),Kp(4), 2.14);
    A.set(K(5),Kp(5),-6.87); 

    ITensor U(Kp),D,V;

    svd(A,U,D,V,{"SVDMethod","gesdd","Truncate",false});

    CHECK(norm(A-U*D*V) < 1E-12);
    }

SECTION("zgesdd test")
    {
    Index m = Index("m",3);
    Index n = Index("n",4);

    ITensor A(m,n);
    A.set(m(1),n(1),-5.40+7.40_i);
    A.set(m(1),n(2), 6.00+6.38_i);
    A.set(m(1),n(3), 9.91+0.16_i);
    A.set(m(1),n(4),-5.28-4.16_i);

    A.set(m(2),n(1), 1.09+1.55_i);
    A.set(m(2),n(2), 2.60+0.07_i);
    A.set(m(2),n(3), 3.98-5.26_i);
    A.set(m(2),n(4), 2.03+1.11_i);

    A.set(m(3),n(1), 9.88+1.91_i);
    A.set(m(3),n(2), 4.92+6.31_i);
    A.set(m(3),n(3),-2.11+7.39_i);
    A.set(m(3),n(4),-9.81-8.98_i);

    ITensor U(m),D,V;

    svd(A,U,D,V,{"SVDMethod","gesdd","Truncate",false});

    CHECK(norm(A-U*D*V) < 1E-12);
    }

SECTION("Transpose dgesdd SVD")
    {
    Index a("a",3),
          b("b",2);

    ITensor A(a,b);

    A.set(a(1),b(1),1);
    A.set(a(2),b(1),2);
    A.set(a(3),b(1),3);
    A.set(a(1),b(2),4);
    A.set(a(2),b(2),5);
    A.set(a(3),b(2),6);

    ITensor U(b),D,V;
    svd(A,U,D,V,{"SVDMethod","gesdd","Truncate",false});

    CHECK(norm(A-U*D*V) < 1E-12);
    }

SECTION("dgesdd SVD")
    {
    Index i("i",2),
          j("j",3),
          k("k",4),
          l("l",5);

    SECTION("Case 1")
        {
        auto T = randomTensor(i,j,k);

        ITensor U(i,j),D,V;

        svd(T,U,D,V,{"SVDMethod","gesdd"});

        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasindex(U,i));
        CHECK(hasindex(U,j));
        CHECK(hasindex(V,k));
        }

    SECTION("Case 2")
        {
        auto T = randomTensor(i,j,k);

        ITensor U(i,k),D,V;
        svd(T,U,D,V,{"SVDMethod","gesdd"});
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasindex(U,i));
        CHECK(hasindex(U,k));
        CHECK(hasindex(V,j));
        }

    SECTION("Case 3")
        {
        auto T = randomTensor(i,k,prime(i));

        ITensor U(i,prime(i)),D,V;
        svd(T,U,D,V,{"SVDMethod","gesdd"});
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasindex(U,i));
        CHECK(hasindex(U,prime(i)));
        CHECK(hasindex(V,k));
        }

    }

SECTION("IQTensor dgesdd SVD")
    {

    SECTION("Regression Test 1")
        {
        //Oct 5, 2015: was encountering a 
        //bad memory access bug with this code
        IQIndex u("u",Index{"u+2",1},QN(+2),
                      Index{"u00",1},QN( 0),
                      Index{"u-2",1},QN(-2));
        IQIndex v("v",Index{"v+2",1},QN(+2),
                      Index{"v00",1},QN( 0),
                      Index{"v-2",1},QN(-2));

        auto S = randomTensor(QN(),u,v);
        IQTensor U(u),D,V;
        svd(S,U,D,V,{"SVDMethod","gesdd"});

        CHECK(norm(S-U*D*V) < 1E-12);
        }

    SECTION("Regression Test 2")
        {
        //Feb 10, 2016: code that fixes sign of
        //singular values to be positive was broken
		auto s1 = IQIndex("s1",Index("s1+",1,Site),QN(+1),Index("s1-",1,Site),QN(-1));
		auto s2 = IQIndex("s2",Index("s2+",1,Site),QN(+1),Index("s2-",1,Site),QN(-1));
		auto sing = IQTensor(s1,s2);
		sing.set(s1(1),s2(2), 1./sqrt(2));
		sing.set(s1(2),s2(1),-1./sqrt(2));
		auto prod = IQTensor(s1,s2);
		prod.set(s1(1),s2(2),1.);
		auto psi = sing*sin(0.1)+prod*cos(0.1);
		psi /= norm(psi);
        psi.scaleTo(-1.);
		IQTensor A(s1),D,B;
		svd(psi,A,D,B,{"SVDMethod","gesdd"});
        CHECK(norm(psi-A*D*B) < 1E-12);
        }

    }

}
