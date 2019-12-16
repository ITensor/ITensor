#include "test.h"
#include "itensor/decomp.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using namespace std;

TEST_CASE("Decomposition Tests")
{

SECTION("SVD interface")
    {
    auto i = Index(3,"i"),
         j = Index(2,"j"),
         k = Index(4,"k");
    auto A = randomITensor(i,j,k);
    auto [U,S,V] = svd(A,i,j);
    CHECK( norm(A-U*S*V) < 1E-12 );
    std::tie(U,S,V) = svd(A,i);
    CHECK( norm(A-U*S*V) < 1E-12 );
    std::tie(U,S,V) = svd(A,i,Args("Cutoff=",1E-13));
    CHECK( norm(A-U*S*V) < 1E-12 );
    }

SECTION("Factor interface")
    {
    auto i = Index(3,"i"),
         j = Index(2,"j"),
         k = Index(4,"k");
    auto A = randomITensor(i,j,k);
    auto [X,Y] = factor(A,i,j);
    CHECK( norm(A-X*Y) < 1E-12 );
    std::tie(X,Y) = factor(A,i);
    CHECK( norm(A-X*Y) < 1E-12 );
    std::tie(X,Y) = factor(A,i,Args("Cutoff=",1E-13));
    CHECK( norm(A-X*Y) < 1E-12 );
    }

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

    auto [U,D,V] = svd(A,{b},{"Truncate",false});

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

    long maxdim = 2*origm,
         mindim = 1;
    Real cutoff = 0;
    //bool absoluteCutoff = false,
    //     doRelCutoff = false;
    Real truncerr = 0,
         docut_lower = 0,
         docut_upper = 0;
    int ndegen = 0;

    SECTION("Case 0")
        {
        //Check that with unrestrictive settings
        //nothing gets truncated
        tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(p,maxdim,mindim,cutoff);
        size_t m = p.size();
        CHECK(m==origm);
        }

    SECTION("Case 1")
        {
        //Check that maxdim is enforced
        maxdim = 10;
        tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(p,maxdim,mindim,cutoff);
        long m = p.size();
        CHECK(m==maxdim);
        }

    SECTION("Case 2")
        {
        //Check that maxdim is enforced
        maxdim = 8;
        tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(p,maxdim,mindim,cutoff);
        long m = p.size();
        CHECK(m==maxdim);
        }

    SECTION("Case 3")
        {
        //Check that mindim is enforced
        mindim = 10;
        cutoff = 0.01;
        tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(p,maxdim,mindim,cutoff);
        long m = p.size();
        CHECK(m==mindim);
        }

    SECTION("Case 4")
        {
        //Check that cutoff is enforced
        //and truncerr is correct
        cutoff = 1E-5;
        auto origp = p;
        tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(p,maxdim,mindim,cutoff);
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
        auto T = randomITensor(i,j,k);

        ITensor U(i,j),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasIndex(U,i));
        CHECK(hasIndex(U,j));
        CHECK(hasIndex(V,k));
        }

    SECTION("Case 2")
        {
        auto T = randomITensor(i,j,k);

        ITensor U(i,k),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasIndex(U,i));
        CHECK(hasIndex(U,k));
        CHECK(hasIndex(V,j));
        }

    SECTION("Case 3")
        {
        auto T = randomITensor(i,k,prime(i));

        ITensor U(i,prime(i)),D,V;
        svd(T,U,D,V);
        CHECK(norm(T-U*D*V) < 1E-12);
        CHECK(hasIndex(U,i));
        CHECK(hasIndex(U,prime(i)));
        CHECK(hasIndex(V,k));
        }

    }

SECTION("QN ITensor SVD")
    {

    SECTION("Regression Test 1")
        {
        //Oct 5, 2015: was encountering a 
        //bad memory access bug with this code
        Index u(QN(+2),1,
                QN( 0),1,
                QN(-2),1);
        Index v(QN(+2),1,
                QN( 0),1,
                QN(-2),1);

        auto S = randomITensor(QN(),u,v);
        ITensor U(u),D,V;
        svd(S,U,D,V);

        CHECK(norm(S-U*D*V) < 1E-12);
        }

    SECTION("Regression Test 2")
        {
        //Feb 10, 2016: code that fixes sign of
        //singular values to be positive was broken
		auto s1 = Index(QN(+1),1,QN(-1),1,"s1");
		auto s2 = Index(QN(+1),1,QN(-1),1,"s2");
		auto sing = ITensor(s1,s2);
		sing.set(s1(1),s2(2), 1./sqrt(2));
		sing.set(s1(2),s2(1),-1./sqrt(2));
		auto prod = ITensor(s1,s2);
		prod.set(s1(1),s2(2),1.);
		auto psi = sing*sin(0.1)+prod*cos(0.1);
		psi /= norm(psi);
        psi.scaleTo(-1.);
		ITensor A(s1),D,B;
		svd(psi,A,D,B);
        CHECK(norm(psi-A*D*B) < 1E-12);
        }

    }

 SECTION("QR Decomposition")
   {
     Index i(3),
       j(4),
       k(5),
       l(6),
       u(8),
       v(5);

     SECTION("Case 1 Complete")
       {
	 auto T = randomITensor(u,v);

	 ITensor Q(u),R;
	 qr(T, Q, R, {"Complete", true, "UpperTriangular", true});
	 CHECK(norm(T-Q*R) < 1E-12);
	 auto qrlink = findIndex(Q.inds(),"Link,QR");
	 CHECK(dim(qrlink) == dim(u));
	 auto QQ = Q*dag(prime(Q,u));
	 auto QQt =  Q*dag(prime(Q,qrlink));
	 for(auto r : range1(dim(u)))
	   for(auto c : range1(dim(u)))
	     {
	       if (r==c)
		 {
		   CHECK_CLOSE(elt(QQ,r,c),1.0);
		   CHECK_CLOSE(elt(QQt,r,c),1.0);
		 }
	       else
		 {
		   CHECK_CLOSE(elt(QQ,r,c),0.0);
		   CHECK_CLOSE(elt(QQt,r,c),0.0);
		   if (c <= dim(v) and r > c)
		     CHECK(elt(R,r,c) == 0); //Check R upper triangular
		 }
	     }	 
       }

     SECTION("Case 1 Thin")
       {
	 auto T = randomITensor(u,v);

	 ITensor Q(u),R;
	 qr(T, Q, R, {"Complete", false, "UpperTriangular", true});
	 CHECK(norm(T-Q*R) < 1E-12);
	 auto qrlink = findIndex(Q.inds(),"Link,QR");
	 CHECK(dim(qrlink) == dim(v));
	 //auto QQ = Q*dag(prime(Q,u)); Thin Q is only right unitary
	 auto QQt =  Q*dag(prime(Q,qrlink));
	 for(auto r : range1(dim(v)))
	   for(auto c : range1(dim(v)))
	     {
	       if (r==c)
		 {
		   //CHECK_CLOSE(elt(QQ,r,c),1.0);
		   CHECK_CLOSE(elt(QQt,r,c),1.0);
		 }
	       else
		 {
		   //CHECK_CLOSE(elt(QQ,r,c),0.0);
		   CHECK_CLOSE(elt(QQt,r,c),0.0);
		   if (r > c)
		     CHECK(elt(R,r,c) == 0); //Check R upper triangular
		 }
	     }	 
       }

     SECTION("Case 1 Not Upper Triangular")
       {
	 auto T = randomITensor(u,v);

	 ITensor Q(u),R;
	 qr(T, Q, R, {"Complete", true, "UpperTriangular", true});
	 CHECK(norm(T-Q*R) < 1E-12);
	 auto qrlink = findIndex(Q.inds(),"Link,QR");
	 CHECK(dim(qrlink) == dim(u));
	 auto QQ = Q*dag(prime(Q,u));
	 auto QQt =  Q*dag(prime(Q,qrlink));
	 for(auto r : range1(dim(u)))
	   for(auto c : range1(dim(u)))
	     {
	       if (r==c)
		 {
		   CHECK_CLOSE(elt(QQ,r,c),1.0);
		   CHECK_CLOSE(elt(QQt,r,c),1.0);
		 }
	       else
		 {
		   CHECK_CLOSE(elt(QQ,r,c),0.0);
		   CHECK_CLOSE(elt(QQt,r,c),0.0);
		 }
	     }
       }

     SECTION("Case 2")
       {
	 auto T = randomITensor(i,j,k);

	 ITensor Q(i,k),R;
	 qr(T, Q, R);
	 CHECK(norm(T-Q*R) < 1E-12);
	 CHECK(hasIndex(Q,i));
	 CHECK(hasIndex(Q,k));
	 CHECK(hasIndex(R,j));
       }

     SECTION("Case 3")
       {
	 auto T = randomITensor(i,k,prime(i));

	 ITensor Q(i,prime(i)),R;
	 qr(T, Q, R);
	 CHECK(norm(T-Q*R) < 1E-12);
	 CHECK(hasIndex(Q,i));
	 CHECK(hasIndex(Q,prime(i)));
	 CHECK(hasIndex(R,k));
       }
     
     SECTION("Case 4")
       {
	 auto T = randomITensor(i,j,k);

	 ITensor Q(i),R;
	 qr(T, Q, R);
	 CHECK(norm(T-Q*R) < 1E-12);
	 CHECK(hasIndex(Q,i));
	 CHECK(hasIndex(R,j));
	 CHECK(hasIndex(R,k));
       }

   }
 
 SECTION("Complex QR Decomposition")
   {
     Index i(3),
       j(4),
       k(5),
       l(6),
       u(8),
       v(5);


     auto T = randomITensorC(u,v);

     ITensor Q(u),R;
     qr(T, Q, R, {"Complete", true, "UpperTriangular", true});
     CHECK(norm(T-Q*R) < 1E-12);
     auto qrlink = findIndex(Q.inds(),"Link,QR");
     CHECK(dim(qrlink) == dim(u));
     auto QQ = Q*dag(prime(Q,u));
     auto QQt =  Q*dag(prime(Q,qrlink));
     for(auto r : range1(dim(u)))
       for(auto c : range1(dim(u)))
	 {
	   if (r==c)
	     {
	       CHECK_CLOSE(eltC(QQ,r,c),1.0);
	       CHECK_CLOSE(eltC(QQt,r,c),1.0);
	     }
	   else
	     {
	       CHECK_CLOSE(eltC(QQ,r,c),0.0);
	       CHECK_CLOSE(eltC(QQt,r,c),0.0);
	       if (c <= dim(v) and r > c)
		 CHECK(eltC(R,r,c) == 0.0); //Check R upper triangular
	     }
	 }
   }
 

 SECTION("QN ITensor QR")
   {
     SECTION("Zero Divergence")
       {
	 Index u(QN(+2),3,
		 QN( 0),2,
		 QN(-2),2);
	 Index v(QN(+2),2,
		 QN( 0),2,
		 QN(-2),1);

	 auto T = randomITensor(QN(),u,v);

	 ITensor Q(u),R;
	 qr(T, Q, R, {"Complete", true, "UpperTriangular", true});
	 CHECK(norm(T-Q*R) < 1E-12);
	 auto qrlink = findIndex(Q.inds(),"Link,QR");
	 CHECK(dim(qrlink) == dim(u));
	 auto QQ = Q*dag(prime(Q,u));
	 auto QQt =  Q*dag(prime(Q,qrlink));
	 for(auto r : range1(dim(u)))
	   for(auto c : range1(dim(u)))
	     {
	       if (r==c)
		 {
		   CHECK_CLOSE(elt(QQ,r,c),1.0);
		   CHECK_CLOSE(elt(QQt,r,c),1.0);
		 }
	       else
		 {
		   CHECK_CLOSE(elt(QQ,r,c),0.0);
		   CHECK_CLOSE(elt(QQt,r,c),0.0);
		   if (c <= dim(v) and r > c)
		     CHECK(elt(R,r,c) == 0); //Check R upper triangular
		 }
	     }
       }
     SECTION("Zero Divergence Thin")
       {
	 Index u(QN(+2),3,
		 QN( 0),2,
		 QN(-2),2);
	 Index v(QN(+2),2,
		 QN( 0),2,
		 QN(-2),1);

	 auto T = randomITensor(QN(), u,v);

	 ITensor Q(u),R;
	 qr(T, Q, R, {"Complete", false, "UpperTriangular", true});
	 CHECK(norm(T-Q*R) < 1E-12);
	 auto qrlink = findIndex(Q.inds(),"Link,QR");
	 CHECK(dim(qrlink) == dim(v));
	 //auto QQ = Q*dag(prime(Q,u)); Thin Q is only right unitary
	 auto QQt =  Q*dag(prime(Q,qrlink));
	 for(auto r : range1(dim(v)))
	   for(auto c : range1(dim(v)))
	     {
	       if (r==c)
		 {
		   //CHECK_CLOSE(elt(QQ,r,c),1.0);
		   CHECK_CLOSE(elt(QQt,r,c),1.0);
		 }
	       else
		 {
		   //CHECK_CLOSE(elt(QQ,r,c),0.0);
		   CHECK_CLOSE(elt(QQt,r,c),0.0);
		   if (r > c)
		     CHECK(elt(R,r,c) == 0); //Check R upper triangular
		 }
	     }	 
       }
     SECTION("Zero Divergence Not Upper Triangular")
       {
	 Index u(QN(+2),3,
		 QN( 0),2,
		 QN(-2),2);
	 Index v(QN(+2),2,
		 QN( 0),2,
		 QN(-2),1);

	 auto T = randomITensor(QN(), u,v);

	 ITensor Q(u),R;
	 qr(T, Q, R, {"Complete", true, "UpperTriangular", true});
	 CHECK(norm(T-Q*R) < 1E-12);
	 auto qrlink = findIndex(Q.inds(),"Link,QR");
	 CHECK(dim(qrlink) == dim(u));
	 auto QQ = Q*dag(prime(Q,u));
	 auto QQt =  Q*dag(prime(Q,qrlink));
	 for(auto r : range1(dim(u)))
	   for(auto c : range1(dim(u)))
	     {
	       if (r==c)
		 {
		   CHECK_CLOSE(elt(QQ,r,c),1.0);
		   CHECK_CLOSE(elt(QQt,r,c),1.0);
		 }
	       else
		 {
		   CHECK_CLOSE(elt(QQ,r,c),0.0);
		   CHECK_CLOSE(elt(QQt,r,c),0.0);
		 }
	     }
       }
     SECTION("Non-zero Divergence")
       {
	 Index u(QN(+2),3,
		 QN( 0),2,
		 QN(-1),2);
	 Index v(QN(+2),2,
		 QN( 0),2,
		 QN(-1),1);

	 auto T = randomITensor(QN(1),u,v);

	 ITensor Q(u),R;
	 qr(T, Q, R, {"Complete", true, "UpperTriangular", true});
	 CHECK(norm(T-Q*R) < 1E-12);
	 auto qrlink = findIndex(Q.inds(),"Link,QR");
	 CHECK(dim(qrlink) == dim(u));
	 auto QQ = Q*dag(prime(Q,u));
	 auto QQt =  Q*dag(prime(Q,qrlink));
	 for(auto r : range1(dim(u)))
	   for(auto c : range1(dim(u)))
	     {
	       if (r==c)
		 {
		   CHECK_CLOSE(elt(QQ,r,c),1.0);
		   CHECK_CLOSE(elt(QQt,r,c),1.0);
		 }
	       else
		 {
		   CHECK_CLOSE(elt(QQ,r,c),0.0);
		   CHECK_CLOSE(elt(QQt,r,c),0.0);
		   if (c <= dim(v) and r > c)
		     CHECK(elt(R,r,c) == 0); //Check R upper triangular
		 }
	     }
       }
   }

SECTION("Polar")
  {
  auto i = Index(2,"i");
  auto j = Index(2,"j");
  auto k = Index(2,"k");
  auto l = Index(2,"l");
  auto A = randomITensor(i,j,k,l);

  auto [U,P] = polar(A,{i,j});

  auto uinds = commonInds(U,P);
  auto [C,c] = combiner(uinds);

  CHECK_CLOSE(norm(U*P-A),0.);

  auto Uc = U*C;

  auto UUdag = Uc*dag(prime(Uc,c));

  for(auto i : range1(dim(c)))
    CHECK_CLOSE(UUdag.elt(i,i),1.);
  }

SECTION("Polar")
  {
  auto i = Index(QN(),2,In,"i");
  auto j = Index(QN(),2,In,"j");
  auto k = Index(QN(),2,Out,"k");
  auto l = Index(QN(),2,Out,"l");
  auto A = randomITensor(QN(),i,j,k,l);

  auto [U,P] = polar(A,{i,j});

  auto uinds = commonInds(U,P);
  auto [C,c] = combiner(uinds);

  CHECK_CLOSE(norm(U*P-A),0.);

  auto Uc = U*C;

  auto UUdag = Uc*dag(prime(Uc,c));

  for(auto i : range1(dim(c)))
    CHECK_CLOSE(UUdag.elt(i,i),1.);
  }

SECTION("QN ITensor denmatDecomp")
    {
    SECTION("Test 1")
        {
        Index S1(QN(+1),1,
                 QN(-1),1);
        Index S2(QN(+1),1,
                 QN(-1),1);
        Index L1(QN(+2),3,
                 QN(+1),4,
                 QN( 0),8,
                 QN(-1),4,
                 QN(-2),2);
        Index L2(QN(+2),4,
                 QN(+1),6,
                 QN( 0),10,
                 QN(-1),4,
                 QN(-2),3);
        Index L3(QN(+2),2,
                 QN( 0),4,
                 QN(-2),2);

        auto A1 = randomITensor(QN(),L1,S1,L2),
             A2 = randomITensor(QN(),dag(L2),S2,L3);

        auto AA = A1*A2;

        AA *= -1./norm(AA);

        auto spec = denmatDecomp(AA,A1,A2,Fromleft,{"Truncate=",true});
 
        auto AAres = A1*A2;

        CHECK(norm(AA-AAres) < 1E-11);

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
        auto T = randomITensor(i,prime(i));
        T += swapTags(T,"0","1");
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
        auto T = randomITensor(i,prime(i),prime(j),j);
        T += swapTags(T,"0","1");
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
        auto T = randomITensorC(i,prime(i));
        T += conj(swapTags(T,"0","1"));
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(norm(T-conj(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 2 - Primes 1 and 2")
        {
        auto i = Index(10);
        auto T = randomITensor(i,prime(i));
        T += swapTags(T,"0","1");
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

        auto T = randomITensor(prime(i),prime(i,2),prime(i,5),prime(i,6));
        T += swapTags(swapTags(T,"1","5"),"2","6");
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(norm(T-U*D*prime(U,4)) < 1E-12);
        }
    }

SECTION("ITensor diagHermitian (with QNs)")
    {
    SECTION("Rank 2")
        {
        auto I = Index(QN(-1),4,QN(+1),4,"I");
        auto T = randomITensor(QN(),dag(I),prime(I));
        T += swapTags(dag(T),"0","1");
        auto [U,D] = diagHermitian(T);
        auto d = commonIndex(U,D);
        CHECK(hasIndex(U,I));
        CHECK(hasInds(D,{d,prime(d)}));
        CHECK(not hasIndex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 2")
        {
        auto I = Index(QN(-1),4,QN(+1),4,"I");
        auto T = randomITensorC(QN(),dag(I),prime(I));
        CHECK(isComplex(T));
        T += dag(swapTags(T,"0","1"));
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,I));
        CHECK(not hasIndex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Complex Rank 4")
        {
        detail::seed_quickran(1);
        auto I = Index(QN(-1),2,
                       QN(+1),2,"I");
        auto J = Index(QN(-2),2,
                       QN(+2),2,"J");
        //auto T = randomITensorC(QN(),dag(I),prime(I),prime(J),dag(J));
        auto T = randomITensorC(QN(),dag(I),dag(J),prime(J),prime(I));
        CHECK(isComplex(T));
        T += dag(swapTags(T,"0","1"));
        T = swapTags(T,"0","1");
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,I));
        CHECK(hasIndex(U,J));
        CHECK(not hasIndex(U,prime(I)));
        CHECK(not hasIndex(U,prime(J)));
        CHECK(norm(T-dag(U)*D*prime(U)) < 1E-12);
        }

    SECTION("Rank 2 - Primes 1 and 2")
        {
        auto I = Index(QN(-1),4,QN(+1),4,"I");
        auto T = randomITensor(QN(),dag(I),prime(I));
        T += dag(swapTags(T,"0","1"));
        //Raise prime level of T
        //and prime level spacing between inds
        T.replaceTags("1","4");
        T.replaceTags("0","1");
        ITensor U,D;
        diagHermitian(T,U,D);
        CHECK(hasIndex(U,prime(I)));
        CHECK(norm(T-dag(U)*D*prime(U,3)) < 1E-12);
        }
    }

SECTION("Truncating (Special Cases)")
  {
  SECTION("All zeros")
    {
    auto vals = vector<Real>({0.0,0.0});
    auto blockdim = vals.size();
    auto nblck = 4;
    auto qns = vector<QNInt>(nblck);
    for( auto b : range(nblck) )
      qns[b] = QNInt(QN(b),blockdim);
    auto iqn = Index(std::move(qns),"i");
    auto iqnp = prime(iqn);
    auto Aqn = ITensor(iqn,dag(iqnp));
    auto blocksize_sum = 0;
    for( auto b : range1(nblck) )
      {
      auto blocksize_b = blocksize(iqn,b);
      for( auto ii : range1(1,blocksize_b) )
        Aqn.set(blocksize_sum+ii,blocksize_sum+ii,vals[ii-1]);
      blocksize_sum += blocksize_b;
      }
    auto [U,S,V] = svd(Aqn,{iqn},{"Cutoff=",1.0});
    auto u = commonIndex(U,S);
    CHECK(dim(u)==1);
    auto v = commonIndex(S,V);
    CHECK(dim(v)==1);
    }

  SECTION("Starting Dimension of One")
    {
    auto iqn = Index(QN(),1,"i");
    auto iqnp = prime(iqn);
    auto Aqn = randomITensor(QN(),iqn,dag(iqnp));
    Aqn /= norm(Aqn);
    auto [U,S,V] = svd(Aqn,{iqn},{"Cutoff=",1.0});
    auto u = commonIndex(U,S);
    CHECK(dim(u)==1);
    auto v = commonIndex(S,V);
    CHECK(dim(v)==1);
    }

  }

SECTION("Truncating Degenerate Values")
  {
  // Make a test Index
  auto v = vector<Real>({2.0,1.0,1.0,1.0,0.5});
  auto blockdim = v.size();
  auto nblck = 4;
  auto qns = vector<QNInt>(nblck);
  for( auto b : range(nblck) )
    qns[b] = QNInt(QN(b),blockdim);
  auto iqn = Index(std::move(qns),"i");

  // Make a test ITensor with degeneracies
  auto iqnp = prime(iqn);
  auto Aqn = ITensor(iqn,dag(iqnp));
  auto blocksize_sum = 0;
  for( auto b : range1(nblck) )
    {
    auto blocksize_b = blocksize(iqn,b);
    for( auto ii : range1(1,blocksize_b) )
      Aqn.set(blocksize_sum+ii,blocksize_sum+ii,v[ii-1]);
    blocksize_sum += blocksize_b;
    }
  Aqn /= norm(Aqn);

  auto i = removeQNs(iqn);
  auto A = removeQNs(Aqn);

  // Value used for truncation
  auto cutoff_dim = 9;

  SECTION("diagPosSemiDef")
    {
    SECTION("MaxDim")
      {
      auto cutoff = 0.0;
      SECTION("No QNs, don't truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(A,{"MaxDim=",cutoff_dim,
                                         "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==cutoff_dim);
        }
      SECTION("No QNs, truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(A,{"RespectDegenerate=",true,
                                         "MaxDim=",cutoff_dim,
                                         "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==nblck);
        }
      SECTION("QNs, don't truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(Aqn,{"MaxDim=",cutoff_dim,
                                           "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==cutoff_dim);
        }
      SECTION("QNs, truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(Aqn,{"RespectDegenerate=",true,
                                           "MaxDim=",cutoff_dim,
                                           "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==nblck);
        }
      }

    SECTION("MinDim")
      {
      auto cutoff = 1.0;
      SECTION("No QNs, don't truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(A,{"MinDim=",cutoff_dim,
                                         "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==cutoff_dim);
        }
      SECTION("No QNs, truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(A,{"RespectDegenerate=",true,
                                         "MinDim=",cutoff_dim,
                                         "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==dim(i)-nblck);
        }
      SECTION("QNs, don't truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(Aqn,{"MinDim=",cutoff_dim,
                                           "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==cutoff_dim);
        }
      SECTION("QNs, truncate degenerate")
        {
        auto [U,D] = diagPosSemiDef(Aqn,{"RespectDegenerate=",true,
                                           "MinDim=",cutoff_dim,
                                           "Cutoff=",cutoff});
        auto u = commonIndex(U,D);
        CHECK(dim(u)==dim(i)-nblck);
        }
      }
    }

  SECTION("svd")
    {
    SECTION("MaxDim")
      {
      auto cutoff = 0.0;
      SECTION("No QNs, don't truncate degenerate")
        {
        auto [U,S,V] = svd(A,{i},{"MaxDim=",cutoff_dim,
                                      "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==cutoff_dim);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==cutoff_dim);
        }
      SECTION("No QNs, truncate degenerate")
        {
        auto [U,S,V] = svd(A,{i},{"RespectDegenerate=",true,
                                      "MaxDim=",cutoff_dim,
                                      "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==nblck);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==nblck);
        }
      SECTION("QNs, don't truncate degenerate")
        {
        auto [U,S,V] = svd(Aqn,{i},{"MaxDim=",cutoff_dim,
                                        "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==cutoff_dim);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==cutoff_dim);
        }
      SECTION("QNs, truncate degenerate")
        {
        auto [U,S,V] = svd(Aqn,{i},{"RespectDegenerate=",true,
                                        "MaxDim=",cutoff_dim,
                                        "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==nblck);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==nblck);
        }
      }

    SECTION("MinDim")
      {
      auto cutoff = 1.0;
      SECTION("No QNs, don't truncate degenerate")
        {
        auto [U,S,V] = svd(A,{i},{"MinDim=",cutoff_dim,
                                      "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==cutoff_dim);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==cutoff_dim);
        }
      SECTION("No QNs, truncate degenerate")
        {
        auto [U,S,V] = svd(A,{i},{"RespectDegenerate=",true,
                                      "MinDim=",cutoff_dim,
                                      "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==dim(i)-nblck);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==dim(i)-nblck);
        }
      SECTION("QNs, don't truncate degenerate")
        {
        auto [U,S,V] = svd(Aqn,{i},{"MinDim=",cutoff_dim,
                                        "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==cutoff_dim);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==cutoff_dim);
        }
      SECTION("QNs, truncate degenerate")
        {
        auto [U,S,V] = svd(Aqn,{i},{"RespectDegenerate=",true,
                                        "MinDim=",cutoff_dim,
                                        "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==dim(i)-nblck);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==dim(i)-nblck);
        }
      }

    SECTION("Cutoff")
      {
      auto cutoff = 0.25;
      auto dim_goal = 10;
      SECTION("No QNs, don't truncate degenerate")
        {
        auto [U,S,V] = svd(A,{i},{"Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==dim_goal);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==dim_goal);
        }
      SECTION("No QNs, truncate degenerate")
        {
        auto [U,S,V] = svd(A,{i},{"RespectDegenerate=",true,
                                      "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==dim(i)-nblck);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==dim(i)-nblck);
        }
      SECTION("QNs, don't truncate degenerate")
        {
        auto [U,S,V] = svd(Aqn,{i},{"Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==dim_goal);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==dim_goal);
        }
      SECTION("QNs, truncate degenerate")
        {
        auto [U,S,V] = svd(Aqn,{i},{"RespectDegenerate=",true,
                                        "Cutoff=",cutoff});
        auto u = commonIndex(U,S);
        CHECK(dim(u)==dim(i)-nblck);
        auto v = commonIndex(S,V);
        CHECK(dim(v)==dim(i)-nblck);
        }
      }

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

    SECTION("QN ITensor case")
        {
        auto s = Index(QN(-1),1,QN(+1),1);
        auto Z = ITensor(dag(s),prime(s));
        Z.set(s(1),prime(s)(1),1);
        Z.set(s(2),prime(s)(2),-1);

        auto Id = ITensor(dag(s),prime(s));
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

SECTION("Eigen")
    {
    auto i = Index(2,"i"),
         j = Index(2,"j"),
         k = Index(2,"k");
    auto I = IndexSet(i,j,k);
    auto Ip = prime(I);

    SECTION("Case 1")
      {
      auto A = randomITensor({Ip,I});

      auto [P,D] = eigen(A,{"Tags=","test"});

      auto d = commonIndex(P,D);

      CHECK(hasTags(d,"test"));
      CHECK_CLOSE(norm(A*P - prime(P)*D),0.);
      }

    SECTION("Case 2")
      {
      auto A = randomITensor({I,Ip});

      auto [P,D] = eigen(A,{"Tags=","test"});

      auto d = commonIndex(P,D);

      CHECK(hasTags(d,"test"));
      CHECK_CLOSE(norm(A*P - prime(P)*D),0.);
      }

    }
}
