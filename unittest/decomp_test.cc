#include "test.h"
#include "itensor/decomp.h"
#include "itensor/util/print_macro.h"

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
        size_t m = p.size();
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

SECTION("ITensor SVD (degeneracy test)")
    {
    Index i("i",3);
    auto T = ITensor(i,prime(i));
    T.set(1,1,1.0);
    T.set(2,2,1.0);
    T.set(3,3,2.0);

    SECTION("Case 1: Ignore degeneracy")
        {
        ITensor U(i),D,V;
        //IgnoreDegeneracy=true is the default
        svd(T,U,D,V,{"Maxm=",2,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==2);
        }

    SECTION("Case 2: Don't ignore degeneracy")
        {
        ITensor U(i),D,V;
        svd(T,U,D,V,{"Maxm=",2,"IgnoreDegeneracy=",false,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==1);
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

SECTION("IQTensor SVD (degeneracy test)")
    {
    auto i = IQIndex("i",Index("i0",3),QN(0),Index("i1",2),QN(1));
    auto A = IQTensor(i,dag(prime(i)));
    A.set(2,2,2.0);
    A.set(4,4,1.0);
    A.set(5,5,1.0);

    SECTION("Case 1: Ignore degeneracy")
        {
        IQTensor U(i),D,V;
        // Degeneracy ignored by default
        svd(A,U,D,V,{"Maxm=",2,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==2);
        }

    SECTION("Case 2: Don't ignore degeneracy")
        {
        IQTensor U(i),D,V;
        svd(A,U,D,V,{"Maxm=",2,"IgnoreDegeneracy=",false,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==1);
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

SECTION("ITensor diagHermitian")
    {
    SECTION("Rank 2")
        {
        auto i = Index("i",10);
        auto T = randomTensor(i,prime(i));
        T += swapPrime(T,0,1);
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasindex(U,i));
        CHECK(not hasindex(U,prime(i)));
        CHECK(norm(T-U*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 4")
        {
        auto i = Index("i",10);
        auto j = Index("i",4);
        auto T = randomTensor(i,prime(i),prime(j),j);
        T += swapPrime(T,0,1);
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasindex(U,i));
        CHECK(hasindex(U,j));
        CHECK(not hasindex(U,prime(i)));
        CHECK(not hasindex(U,prime(j)));
        CHECK(norm(T-U*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 2")
        {
        auto i = Index("i",10);
        auto T = randomTensorC(i,prime(i));
        T += conj(swapPrime(T,0,1));
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(norm(T-conj(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 2 - Primes 1 and 2")
        {
        auto i = Index("i",10);
        auto T = randomTensor(i,prime(i));
        T += swapPrime(T,0,1);
        //Raise prime level of T
        T.prime();
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasindex(U,prime(i)));
        CHECK(norm(T-U*D*prime(U)) < 1E-12);
        }

    SECTION("Multiple Prime Levels")
        {
        auto i = Index("i",3);

        auto T = randomTensor(prime(i),prime(i,2),prime(i,5),prime(i,6));
        T += swapPrime(swapPrime(T,1,5),2,6);
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(norm(T-U*D*prime(U,4)) < 1E-12);
        }

    SECTION("Ignore degeneracy")
        {
        Index i("i",3);
        auto T = ITensor(i,prime(i));
        T.set(1,1,1.0);
        T.set(2,2,1.0);
        T.set(3,3,2.0);
        ITensor U(i),D;
        //IgnoreDegeneracy=true is the default
        diagHermitian(T,U,D,{"Maxm=",2,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==2);
        }

    SECTION("Don't ignore degeneracy")
        {
        Index i("i",3);
        auto T = ITensor(i,prime(i));
        T.set(1,1,1.0);
        T.set(2,2,1.0);
        T.set(3,3,2.0);
        ITensor U(i),D;
        diagHermitian(T,U,D,{"Maxm=",2,"IgnoreDegeneracy=",false,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==1);
        }

    }

SECTION("IQTensor diagHermitian")
    {
    SECTION("Rank 2")
        {
        auto I = IQIndex("I",Index("i-",4),QN(-1),Index("i+",4),QN(+1));
        auto T = randomTensor(QN(),dag(I),prime(I));
        T += dag(swapPrime(T,0,1));
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasindex(U,I));
        CHECK(not hasindex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 2")
        {
        auto I = IQIndex("I",Index("i-",4),QN(-1),Index("i+",4),QN(+1));
        auto T = randomTensorC(QN(),dag(I),prime(I));
        CHECK(isComplex(T));
        T += dag(swapPrime(T,0,1));
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasindex(U,I));
        CHECK(not hasindex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 4")
        {
        detail::seed_quickran(1);
        auto I = IQIndex("I",Index("i-",2),QN(-1),
                             Index("i+",2),QN(+1));
        auto J = IQIndex("J",Index("j-2",2),QN(-2),
                             Index("j+2",2),QN(+2));
        //auto T = randomTensorC(QN(),dag(I),prime(I),prime(J),dag(J));
        auto T = randomTensorC(QN(),dag(I),dag(J),prime(J),prime(I));
        CHECK(isComplex(T));
        T += dag(swapPrime(T,0,1));
        T = swapPrime(T,0,1);
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasindex(U,I));
        CHECK(hasindex(U,J));
        CHECK(not hasindex(U,prime(I)));
        CHECK(not hasindex(U,prime(J)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 2 - Primes 1 and 2")
        {
        auto I = IQIndex("I",Index("i-",4),QN(-1),Index("i+",4),QN(+1));
        auto T = randomTensor(QN(),dag(I),prime(I));
        T += dag(swapPrime(T,0,1));
        //Raise prime level of T
        //and prime level spacing between inds
        T.mapprime(1,4);
        T.mapprime(0,1);
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasindex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U,3)) < 1E-12);
        }

    SECTION("Ignore degeneracy")
        {
        auto i = IQIndex("i",Index("i0",3),QN(0),Index("i1",2),QN(1));
        auto A = IQTensor(i,dag(prime(i)));
        A.set(2,2,2.0);
        A.set(4,4,1.0);
        A.set(5,5,1.0);
        IQTensor U(i),D;
        // Degeneracy ignored by default
        diagHermitian(A,U,D,{"Maxm=",2,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==2);
        }

    SECTION("Don't ignore degeneracy")
        {
        auto i = IQIndex("i",Index("i0",3),QN(0),Index("i1",2),QN(1));
        auto A = IQTensor(i,dag(prime(i)));
        A.set(2,2,2.0);
        A.set(4,4,1.0);
        A.set(5,5,1.0);
        IQTensor U(i),D;
        diagHermitian(A,U,D,{"Maxm=",2,"IgnoreDegeneracy=",false,"Cutoff=",0.0});
        auto l = commonIndex(U,D);
        CHECK(l.m()==1);
        }

    }

SECTION("Exp Hermitian")
    {
    SECTION("ITensor case")
        {
        auto s = Index("s",2);
        auto X = ITensor(s,prime(s));
        X.set(s(1),prime(s)(2),1);
        X.set(s(2),prime(s)(1),1);

        auto Id = ITensor(s,prime(s));
        Id.set(s(1),prime(s)(1),1);
        Id.set(s(2),prime(s)(2),1);

        SECTION("Real tau")
            {
            Real tau = -0.2342;
            auto expX = expHermitian(X,tau);
            CHECK(norm(expX - (cosh(tau)*Id+sinh(tau)*X)) < 1E-12);
            CHECK(not isComplex(expX));
            }

        SECTION("Imag tau")
            {
            Cplx tau = -0.2342_i;
            auto expX = expHermitian(X,tau);
            CHECK(norm(expX - (cos(tau.imag())*Id+1_i*sin(tau.imag())*X)) < 1E-12);
            CHECK(isComplex(expX));
            }
        }

    SECTION("IQTensor case")
        {
        auto s = IQIndex("S",Index("s-",1),QN(-1),
                             Index("s+",1),QN(+1));
        auto Z = IQTensor(dag(s),prime(s));
        Z.set(s(1),prime(s)(1),1);
        Z.set(s(2),prime(s)(2),-1);

        auto Id = IQTensor(dag(s),prime(s));
        Id.set(s(1),prime(s)(1),1);
        Id.set(s(2),prime(s)(2),1);

        SECTION("Real tau")
            {
            Real tau = -0.2342;
            auto expZ = expHermitian(Z,tau);
            CHECK(norm(expZ - (cosh(tau)*Id+sinh(tau)*Z)) < 1E-12);
            CHECK(not isComplex(expZ));
            }

        SECTION("Imag tau")
            {
            Cplx tau = -0.2342_i;
            auto expZ = expHermitian(Z,tau);
            CHECK(norm(expZ - (cos(tau.imag())*Id+1_i*sin(tau.imag())*Z)) < 1E-12);
            CHECK(isComplex(expZ));
            }
        }
    }

}
