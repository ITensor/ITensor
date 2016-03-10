#include "test.h"
#include "itensor/decomp.h"

using namespace itensor;
using namespace std;

TEST_CASE("Decomposition Tests")
{

SECTION("Transpose SVD")
    {
    Index a("a",3),
          b("b",2);

    ITensor A(a,b);
    auto el = 1./sqrt(2);

    A.set(a(1),b(1),el);
    A.set(a(2),b(1),-el);
    A.set(a(3),b(1),el);
    A.set(a(1),b(2),el);
    A.set(a(2),b(2),el);
    A.set(a(3),b(2),el);

    ITensor U(b),D,V;
    svd(A,U,D,V,{"Truncate",false});

    CHECK(norm(A-U*D*V) < 1E-12);
    }

SECTION("Truncate Test")
    {
    size_t origm = 20;
    auto p = Vector(origm);
    
    for(auto n : range(p))
        {
        p(n) = exp(-2.*n)*(1+0.1*detail::quickran());
        }
    p(1) = p(2);
    p(3) = p(5);
    p(4) = p(5);
    p /= sumels(p);
    //Print(p);

    long maxm = 2*origm,
         minm = 1;
    Real cutoff = 0;
    //bool absoluteCutoff = false,
    //     doRelCutoff = false;
    Real truncerr = 0,
         docut = 0;

    SECTION("Case 0")
        {
        //Check that with unrestrictive settings
        //nothing gets truncated
        tie(truncerr,docut) = truncate(p,maxm,minm,cutoff);
        long m = p.size();
        CHECK(m==origm);
        }

    SECTION("Case 1")
        {
        //Check that maxm is enforced
        maxm = 10;
        tie(truncerr,docut) = truncate(p,maxm,minm,cutoff);
        long m = p.size();
        CHECK(m==maxm);
        }

    SECTION("Case 2")
        {
        //Check that maxm is enforced
        maxm = 8;
        tie(truncerr,docut) = truncate(p,maxm,minm,cutoff);
        long m = p.size();
        CHECK(m==maxm);
        }

    SECTION("Case 3")
        {
        //Check that minm is enforced
        minm = 10;
        cutoff = 0.01;
        tie(truncerr,docut) = truncate(p,maxm,minm,cutoff);
        long m = p.size();
        CHECK(m==minm);
        }

    SECTION("Case 4")
        {
        //Check that cutoff is enforced
        //and truncerr is correct
        cutoff = 1E-5;
        auto origp = p;
        tie(truncerr,docut) = truncate(p,maxm,minm,cutoff);
        long m = p.size();
        Real te_check = 0.;
        for(auto n = m; n < long(origp.size()); ++n)
            {
            te_check += origp(n);
            }
        //printfln("truncerr = %.5E",truncerr);
        //printfln("te_check = %.5E",te_check);
        CHECK_CLOSE(truncerr,te_check);
        CHECK(truncerr < cutoff);
        }
    }

SECTION("ITensor SVD")
    {
    Index i("i",3),
          j("j",4),
          k("k",5),
          l("l",6);

    SECTION("Case 1")
        {
        auto T = randomTensor(i,j,k);

        ITensor U(i,j),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasindex(U,i));
        CHECK(hasindex(U,j));
        CHECK(hasindex(V,k));
        }

    SECTION("Case 2")
        {
        auto T = randomTensor(i,j,k);

        ITensor U(i,k),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasindex(U,i));
        CHECK(hasindex(U,k));
        CHECK(hasindex(V,j));
        }

    SECTION("Case 3")
        {
        auto T = randomTensor(i,k,prime(i));

        ITensor U(i,prime(i)),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasindex(U,i));
        CHECK(hasindex(U,prime(i)));
        CHECK(hasindex(V,k));
        }

    }

SECTION("IQTensor SVD")
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
        svd(S,U,D,V);

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
		svd(psi,A,D,B);
        CHECK(norm(psi-A*D*B) < 1E-12);
        }

    }

SECTION("IQTensor denmatDecomp")
    {
    SECTION("Test 1")
        {
        IQIndex S1("S1",Index("s1+",1,Site),QN(+1),
                        Index("s1-",1,Site),QN(-1));
        IQIndex S2("S2",Index("s2+",1,Site),QN(+1),
                        Index("s2-",1,Site),QN(-1));
        IQIndex L1("L1",Index("l1+2",3),QN(+2),
                        Index("l1+1",4),QN(+1),
                        Index("l1 0",8),QN( 0),
                        Index("l1-1",4),QN(-1),
                        Index("l1-2",2),QN(-2));
        IQIndex L2("L2",Index("l2+2",4),QN(+2),
                        Index("l2+1",6),QN(+1),
                        Index("l2 0",10),QN( 0),
                        Index("l2-1",4),QN(-1),
                        Index("l2-2",3),QN(-2));
        IQIndex L3("L3",Index("l3+2",2),QN(+2),
                        Index("l3 0",4),QN( 0),
                        Index("l3-2",2),QN(-2));

        auto A1 = randomTensor(QN(),L1,S1,L2),
             A2 = randomTensor(QN(),dag(L2),S2,L3);

        auto AA = A1*A2;
        AA *= -1./norm(AA);
        auto spec = denmatDecomp(AA,A1,A2,Fromleft);

        CHECK(norm(AA-A1*A2) < 1E-11);

        for(auto eig : spec.eigsKept())
            {
            CHECK(eig >= 0.);
            }
        }
    }

}
