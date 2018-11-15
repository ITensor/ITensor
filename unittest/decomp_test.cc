#include "test.h"
#include "itensor/decomp.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using namespace std;

TEST_CASE("Decomposition Tests")
{

SECTION("Transpose SVD")
    {
    Index a(3),
          b(2);

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
    Index i(3),
          j(4),
          k(5),
          l(6);

    SECTION("Case 1")
        {
        auto T = randomTensor(i,j,k);

        ITensor U(i,j),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasIndex(U,i));
        CHECK(hasIndex(U,j));
        CHECK(hasIndex(V,k));
        }

    SECTION("Case 2")
        {
        auto T = randomTensor(i,j,k);

        ITensor U(i,k),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasIndex(U,i));
        CHECK(hasIndex(U,k));
        CHECK(hasIndex(V,j));
        }

    SECTION("Case 3")
        {
        auto T = randomTensor(i,k,prime(i));

        ITensor U(i,prime(i)),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasIndex(U,i));
        CHECK(hasIndex(U,prime(i)));
        CHECK(hasIndex(V,k));
        }

    }

SECTION("IQTensor SVD")
    {

    SECTION("Regression Test 1")
        {
        //Oct 5, 2015: was encountering a 
        //bad memory access bug with this code
        IQIndex u(Index{1},QN(+2),
                  Index{1},QN( 0),
                  Index{1},QN(-2));
        IQIndex v(Index{1},QN(+2),
                  Index{1},QN( 0),
                  Index{1},QN(-2));

        auto S = randomTensor(QN(),u,v);
        IQTensor U(u),D,V;
        svd(S,U,D,V);

        CHECK(norm(S-U*D*V) < 1E-12);
        }

    SECTION("Regression Test 2")
        {
        //Feb 10, 2016: code that fixes sign of
        //singular values to be positive was broken
		auto s1 = IQIndex(Index(1),QN(+1),Index(1),QN(-1));
		auto s2 = IQIndex(Index(1),QN(+1),Index(1),QN(-1));
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
        IQIndex S1(Index(1),QN(+1),
                   Index(1),QN(-1));
        IQIndex S2(Index(1),QN(+1),
                   Index(1),QN(-1));
        IQIndex L1(Index(3),QN(+2),
                   Index(4),QN(+1),
                   Index(8),QN( 0),
                   Index(4),QN(-1),
                   Index(2),QN(-2));
        IQIndex L2(Index(4),QN(+2),
                   Index(6),QN(+1),
                   Index(10),QN( 0),
                   Index(4),QN(-1),
                   Index(3),QN(-2));
        IQIndex L3(Index(2),QN(+2),
                   Index(4),QN( 0),
                   Index(2),QN(-2));

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
        auto i = Index(10);
        auto T = randomTensor(i,prime(i));
        T += swapPrime(T,0,1);
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,i));
        CHECK(not hasIndex(U,prime(i)));
        CHECK(norm(T-U*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 4")
        {
        auto i = Index(10);
        auto j = Index(4);
        auto T = randomTensor(i,prime(i),prime(j),j);
        T += swapPrime(T,0,1);
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,i));
        CHECK(hasIndex(U,j));
        CHECK(not hasIndex(U,prime(i)));
        CHECK(not hasIndex(U,prime(j)));
        CHECK(norm(T-U*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 2")
        {
        auto i = Index(10);
        auto T = randomTensorC(i,prime(i));
        T += conj(swapPrime(T,0,1));
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(norm(T-conj(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 2 - Primes 1 and 2")
        {
        auto i = Index(10);
        auto T = randomTensor(i,prime(i));
        T += swapPrime(T,0,1);
        //Raise prime level of T
        T.prime();
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,prime(i)));
        CHECK(norm(T-U*D*prime(U)) < 1E-12);
        }

    SECTION("Multiple Prime Levels")
        {
        auto i = Index(3);

        auto T = randomTensor(prime(i),prime(i,2),prime(i,5),prime(i,6));
        T += swapPrime(swapPrime(T,1,5),2,6);
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(norm(T-U*D*prime(U,4)) < 1E-12);
        }
    }

SECTION("IQTensor diagHermitian")
    {
    SECTION("Rank 2")
        {
        auto I = IQIndex(Index(4),QN(-1),Index(4),QN(+1));
        auto T = randomTensor(QN(),dag(I),prime(I));
        T += dag(swapPrime(T,0,1));
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,I));
        CHECK(not hasIndex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 2")
        {
        auto I = IQIndex(Index(4),QN(-1),Index(4),QN(+1));
        auto T = randomTensorC(QN(),dag(I),prime(I));
        CHECK(isComplex(T));
        T += dag(swapPrime(T,0,1));
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,I));
        CHECK(not hasIndex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 4")
        {
        detail::seed_quickran(1);
        auto I = IQIndex(Index(2),QN(-1),
                         Index(2),QN(+1));
        auto J = IQIndex(Index(2),QN(-2),
                         Index(2),QN(+2));
        //auto T = randomTensorC(QN(),dag(I),prime(I),prime(J),dag(J));
        auto T = randomTensorC(QN(),dag(I),dag(J),prime(J),prime(I));
        CHECK(isComplex(T));
        T += dag(swapPrime(T,0,1));
        T = swapPrime(T,0,1);
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,I));
        CHECK(hasIndex(U,J));
        CHECK(not hasIndex(U,prime(I)));
        CHECK(not hasIndex(U,prime(J)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 2 - Primes 1 and 2")
        {
        auto I = IQIndex(Index(4),QN(-1),Index(4),QN(+1));
        auto T = randomTensor(QN(),dag(I),prime(I));
        T += dag(swapPrime(T,0,1));
        //Raise prime level of T
        //and prime level spacing between inds
        T.mapPrime(1,4);
        T.mapPrime(0,1);
        IQTensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U,3)) < 1E-12);
        }
    }

SECTION("Exp Hermitian")
    {
    SECTION("ITensor case")
        {
        auto s = Index(2);
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
        auto s = IQIndex(Index(1),QN(-1),
                         Index(1),QN(+1));
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
