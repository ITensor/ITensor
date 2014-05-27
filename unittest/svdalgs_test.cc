#include "test.h"
#include "svdalgs.h"

using namespace itensor;
using namespace std;

TEST_CASE("SVDAlgsTest")
{
IQIndex S1,S2,L1,L2,Mid,Mid10;


IQTensor Phi0,L,R,Psi;
IQTensor V;
ITensor phi0,l,r,psi;
ITensor v;

Index s1u("Site1 Up",1,Site);
Index s1d("Site1 Dn",1,Site);
Index s2u("Site2 Up",1,Site);
Index s2d("Site2 Dn",1,Site);
Index l1u("Link1 Up",10,Link);
Index l10("Link1 Z0",20,Link);
Index l1d("Link1 Dn",10,Link);
Index l2u("Link2 Up",5,Link);
Index l20("Link2 Z0",50,Link);
Index l2d("Link2 Dn",5,Link);
Index mid("mid");
Index mid10u("mid10u",2);
Index mid10z("mid10z",5);
Index mid10d("mid10d",3);

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
             l2u,QN(+1),
             l20,QN( 0),
             l2d,QN(-1),
             Out);

Mid = IQIndex("Mid",mid,QN(),Out);

Mid10 = IQIndex("Mid10",
                mid10u,QN(+1),
                mid10z,QN( 0),
                mid10d,QN(-1),
                Out);

//Construct a divergenceless IQTensor phi0
Phi0 = IQTensor(L1,S1,S2,L2);
    {
    ITensor uudd(l1u,s1u,s2d,l2d);
    uudd.randomize();
    Phi0 += uudd;

    ITensor zudz(l10,s1u,s2d,l20);
    zudz.randomize();
    Phi0 += zudz;

    ITensor duud(l1d,s1u,s2u,l2d);
    duud.randomize();
    Phi0 += duud;
    }
Phi0 *= 1.0/Phi0.norm();

const QN Zero;
CHECK_EQUAL(div(Phi0),Zero);

phi0 = Phi0.toITensor();

L = IQTensor(L1,S1,conj(Mid10));
    {
    ITensor zuu(l10,s1u,mid10u);
    zuu.randomize();
    L += zuu;

    ITensor udz(l1u,s1d,mid10z);
    udz.randomize();
    L += udz;

    ITensor duz(l1d,s1u,mid10z);
    duz.randomize();
    L += duz;
    }
L *= 1./L.norm();
l = L.toITensor();

R = IQTensor(Mid10,S2,L2);
    {
    ITensor zud(mid10z,s2u,l2d);
    zud.randomize();
    R += zud;

    ITensor udz(mid10u,s2d,l20);
    udz.randomize();
    R += udz;

    ITensor duz(mid10d,s2u,l20);
    duz.randomize();
    R += duz;

    ITensor zdu(mid10z,s2d,l2u);
    zdu.randomize();
    R += zdu;
    }
R *= 1./R.norm();
r = R.toITensor();

Psi = L;
Psi *= R;

psi = Psi.toITensor();

SECTION("SiteSVD")
    {
    //
    //ITensor version
    //
    ITensor a(L1,S1),b(S2,L2);

    //svd.showeigs(true);
    phi0 *= -1.2324;
    Spectrum spec = denmatDecomp(phi0,a,b,Fromleft);

    //Print(((a*b)-phi0).norm());
    CHECK(((a*b)-phi0).norm() < 1E-12 );

    CHECK(spec.truncerr() < 1E-12);

    //
    //IQTensor version
    //

    IQTensor A(L1,S1),B(S2,L2);

    Phi0 *= -10.23432;
    spec = denmatDecomp(Phi0,A,B,Fromleft);

    //Print(((A*B)-Phi0).norm());
    CHECK(((A*B)-Phi0).norm() < 1E-12 );

    CHECK(spec.truncerr() < 1E-12);
    }

SECTION("BondSVD")
    {
    //
    //ITensor version
    //

    ITensor l(L1,S1,Mid),r(Mid,S2,L2);
    ITensor v(Mid);

    phi0 *= -0.235;

    Spectrum spec = csvd(phi0,l,v,r);

    //Print(((l*v*r)-phi0).norm());
    CHECK(((l*v*r)-phi0).norm() < 1E-12 );


    CHECK(spec.truncerr() < 1E-12);

    //
    //IQTensor version
    //

    IQTensor L(L1,S1,Mid),R(Mid,S2,L2);
    IQTensor V(Mid);

    spec = csvd(Phi0,L,V,R);
    
    IQTensor nPhi = L*V*R;

    CHECK((nPhi-Phi0).norm() < 1E-12 );

    CHECK(spec.truncerr() < 1E-12);
    
    }

SECTION("CSVDNorm")
    {
    //
    //ITensor version
    //

    Index rr("rr",1,Link),s1("s1",3,Site),s2("s2",3,Site);
    phi0 = ITensor(s1,s2,rr);

    ITensor l(s1,Mid),r(Mid,s2,rr);
    ITensor v(Mid);

    phi0(s1(3),s2(1),rr(1)) = -0.172148;
    phi0(s1(2),s2(2),rr(1)) = 0.427132;
    phi0(s1(1),s2(3),rr(1)) = -0.88765;
    phi0.scaleTo(-2.127145);

    svd(phi0,l,v,r);
    Index nmr = commonIndex(r,v,Link);
    //ITensor rhoRsvd1 = r * conj(primeind(r,rr));
    //ITensor rhoRsvd2 = r * conj(primeind(r,nmr));
    //PrintDat(r);
    //PrintDat(rhoRsvd1);
    //PrintDat(rhoRsvd2);


    csvd(phi0,l,v,r);

    CHECK_CLOSE(1,l.norm(),1E-2);
    CHECK_CLOSE(1,r.norm(),1E-2);
    CHECK((l*v*r-phi0).norm() < 1E-10);

    //cerr << format("l.norm() = %.3E\n") % l.norm();
    //cerr << format("v.norm() = %.3E\n") % v.norm();
    //cerr << format("r.norm() = %.3E\n") % r.norm();

    //cerr << format("(l*v*r-phi0).norm() = %.3E\n") % (l*v*r-phi0).norm();

    //PrintDat(v * r);
    //ITensor ur = v*r; 
    //ITensor rhoR = ur * conj(primeind(ur,rr));
    //PrintDat(rhoR);

    }

SECTION("AbsoluteCutoff")
    {
    OptSet opts;
    opts.add("AbsoluteCutoff",true);

    //
    //ITensor version
    //

    ITensor a(L1,S1,Mid),b(Mid,S2,L2);
    ITensor c(Mid);

    Real cutoff = 1E-3;
    Spectrum spec = csvd(phi0,a,c,b,opts & Opt("Cutoff",cutoff));
    CHECK(spec.eigsKept()(spec.numEigsKept()) > cutoff);

    cutoff = 1E-5;
    spec = csvd(phi0,a,c,b,opts & Opt("Cutoff",cutoff));
    CHECK(spec.eigsKept()(spec.numEigsKept()) > cutoff);

    cutoff = 1E-7;
    spec = csvd(phi0,a,c,b,opts & Opt("Cutoff",cutoff));
    CHECK(spec.eigsKept()(spec.numEigsKept()) > cutoff);

    //
    //IQTensor version
    //

    IQTensor A(L1,S1,Mid),B(Mid,S2,L2);
    IQTensor C(Mid);

    cutoff = 1E-3;
    spec = csvd(Phi0,A,C,B,opts & Opt("Cutoff",cutoff));
    CHECK(spec.eigsKept()(spec.numEigsKept()) > cutoff);

    cutoff = 1E-5;
    spec = csvd(Phi0,A,C,B,opts & Opt("Cutoff",cutoff));
    CHECK(spec.eigsKept()(spec.numEigsKept()) > cutoff);

    cutoff = 1E-7;
    spec = csvd(Phi0,A,C,B,opts & Opt("Cutoff",cutoff));
    CHECK(spec.eigsKept()(spec.numEigsKept()) > cutoff);
    }

/*
SECTION("UseOrigM")
    {
    SVDWorker svd;
    svd.cutoff(1E-2);
    svd.maxm(6);
    svd.useOrigM(true);
    //svd.showeigs(true);

    //
    //IQTensor version
    //

    svd.csvd(Psi,L,V,R);

    CHECK_EQUAL(Mid10.m(),V.index(1).m());

    IQTensor nPsi = L*V*R;

    CHECK((nPsi-Psi).norm() < 1E-12);
    CHECK(svd.truncerr() < 1E-12);

    //
    //ITensor version
    //

    svd.csvd(psi,l,v,r);

    CHECK_EQUAL(Mid10.m(),v.index(1).m());
    
    CHECK(((l*v*r)-psi).norm() < 1E-12 );
    CHECK(svd.truncerr() < 1E-12);
    }
    */

/*
SECTION("SvdRank2")
    {
    int n = 10, m = 5;

    Matrix M(n,m);
    Matrix dd(n,n); dd = 0.0;
    for(int i = 1; i <= n; i++)
        {
        Real eig = pow(0.5,5*(i-1));
        dd(i,i) = eig;
        }
    Matrix uu(n,n), vv(n,m);
    uu.randomize(); vv.randomize();
    Orthog(uu,n,2);
    Orthog(vv,n,2);
    M = uu * dd * vv;

    //Matrix UU,VV;
    //Vector DD;
    //SVD(M,UU,DD,VV);
    //Print(DD);

    Index ui("ui",n), vi("vi",m);
    ITensor T(ui,vi,M);
    ITensor U(ui), V(vi);

    SVDWorker svd;
    svd.showeigs(true);

    ITensor D;
    svd.svdRank2(T,U,D,V);
    ITensor nT = U * D * V;

    cerr << format("\n(T-nT).norm() = %.3E\n") % ((T-nT).norm()) << endl;

    //
    //IQTensor version
    //

    IQIndex uI("uI",ui,QN(0),Out), vI("vI",vi,QN(0),Out);
    IQTensor TT(uI,vI);
    TT += T;

    IQTensor UU(uI), VV(vI);
    IQTensor DD;

    svd.svdRank2(TT,UU,DD,VV);

    }
    */

SECTION("ComplexSVD")
    {
    IQTensor TR(L1(1),prime(L1(31))),
             TI(L1(1),prime(L1(31)));
    IQTensor T = Complex_1*TR+Complex_i*TI;
    T.randomize();
    T *= 1.0/T.norm();

    IQTensor U(L1),D,V;
    svd(T,U,D,V);

    IQTensor diff = T-U*D*V;

    CHECK(diff.norm() < 1E-10);
    }

SECTION("ComplexDenmat")
    {
    Index r("r",4),c("c",4);
    ITensor rr(r,c),ri(r,c);

    rr(r(1),c(1)) = 0.2;
    rr(r(2),c(2)) = 0.2;
    rr(r(2),c(3)) = -0.00492164;
    rr(r(3),c(2)) = -0.00492164;
    rr(r(3),c(3)) = 0.2;
    rr(r(4),c(4)) = 0.4;

    ri(r(2),c(3)) = 0.0983495255;
    ri(r(3),c(2)) = -0.0983495255; 

    ITensor rho = rr*Complex_1+ri*Complex_i;
    rho.scaleTo(1);
    //PrintDat(rho);

    ITensor U(r),V;
    ITensor D;

    svd(rho,U,D,V);
    //PrintDat(U);
    //PrintDat(D);
    //PrintDat(V);
    }

SECTION("Diagonalization")
    {
    Index s1("s1",2,Site),s2("s2",2,Site);
    ITensor M(s1,s2,prime(s2),prime(s1));
    M.randomize();
    M = M + swapPrime(M,0,1);
    M *= 0.5;

    ITensor U;
    ITensor D;
    diagHermitian(M,U,D);

    CHECK((M-(prime(U)*D*conj(U))).norm() < 1E-14);

    //////////////////////////

    IQTensor T(conj(S2)(1),conj(S1)(2),prime(S1)(1),prime(S2)(2));
    T.randomize();
    T = T + swapPrime(T,0,1);
    T *= 0.5;
    //IQTensor T(conj(S1),prime(S1));
    //T(conj(S1)(1),prime(S1)(1)) = 1;
    //T(conj(S1)(2),prime(S1)(2)) = -1;

    IQTensor UU;
    IQTensor DD;
    diagHermitian(T,UU,DD);

    CHECK((T-(prime(UU)*DD*conj(UU))).norm() < 1E-14);
    }

SECTION("ComplexDiagonalization")
    {
    Index s1("s1",2,Site),s2("s2",2,Site);
    ITensor Mr(s1,s2,prime(s2),prime(s1)),
            Mi(s1,s2,prime(s2),prime(s1));
    Mr.randomize();
    Mi.randomize();
    ITensor M = Complex_1*Mr + Complex_i*Mi;
    M = M + conj(swapPrime(M,0,1));
    M *= 0.5;

    ITensor U;
    ITensor D;
    diagHermitian(M,U,D);

    CHECK((M-(prime(U)*D*conj(U))).norm() < 1E-14);

    //////////////////////////

    IQTensor Tr(conj(S2)(1),conj(S1)(2),prime(S1)(1),prime(S2)(2)),
             Ti(conj(S2)(1),conj(S1)(2),prime(S1)(1),prime(S2)(2));
    Tr.randomize();
    Ti.randomize();
    IQTensor T = Complex_1*Tr + Complex_i*Ti;
    T = T + conj(swapPrime(T,0,1));
    T *= 0.5;

    IQTensor UU;
    IQTensor DD;
    diagHermitian(T,UU,DD);

    CHECK((T-(prime(UU)*DD*conj(UU))).norm() < 1E-14);

    //Arrows the other way
    diagHermitian(conj(T),UU,DD);

    CHECK((conj(T)-(prime(UU)*DD*conj(UU))).norm() < 1E-14);
    }

SECTION("OrthoDecomp")
    {
    ITensor phi(Phi0);
    ITensor cphi(phi);
    phi.randomize();
    cphi.randomize();
    phi += Complex_i*cphi;

    ITensor A(L1,S1),B;
    orthoDecomp(phi,A,B,Fromleft);

    CHECK((phi-A*B).norm() < 1E-12);

    //Check that A is orthogonal
    CHECK(!A.isComplex());
    Index mid = commonIndex(A,B);
    ITensor Id(prime(mid),mid,1);
    CHECK((Id-A*prime(A,mid)).norm() < 1E-12);


    //Other direction

    orthoDecomp(phi,A,B,Fromright);

    CHECK((phi-A*B).norm() < 1E-12);

    //Check that B is orthogonal
    CHECK(!B.isComplex());
    mid = commonIndex(A,B);
    Id = ITensor(prime(mid),mid,1);
    CHECK((Id-B*prime(B,mid)).norm() < 1E-12);


    //
    // IQTensor version
    //

    Phi0.randomize();
    IQTensor cPhi0(Phi0);
    cPhi0.randomize();
    Phi0 += Complex_i*cPhi0;
    IQTensor C(L1,S1),D;
    orthoDecomp(Phi0,C,D,Fromleft);

    CHECK((Phi0-C*D).norm() < 1E-12);

    CHECK(!C.isComplex());
    }


SECTION("EigDecomp")
    {
    Index i("i",4),
          j("j",10);

    ITensor Mr(i,j,prime(i),prime(j)),
            Mi(i,j,prime(i),prime(j));

    Mr.randomize();
    Mi.randomize();

    ITensor V,D;
    eigDecomp(Mr,V,D);

    CHECK((Mr*V - prime(V)*D).norm() < 1E-11);

    ITensor M = Mr+Complex_i*Mi;
    eigDecomp(M,V,D);

    CHECK((M*V - prime(V)*D).norm() < 1E-11);
    }

SECTION("IQEigDecomp")
    {
    IQIndex I("J",Index("I-",2),QN(-1),
                  Index("I0",2),QN(+0),
                  Index("I+",2),QN(+1));
    IQIndex J("J",Index("J-",2),QN(-1),
                  Index("J0",2),QN(+0),
                  Index("J+",2),QN(+1));

    IQTensor Mr(conj(I)(1),conj(J)(J.m()),prime(I(1)),prime(J(J.m()))),
             Mi(conj(I)(1),conj(J)(J.m()),prime(I(1)),prime(J(J.m())));

    Mr.randomize();
    Mi.randomize();

    IQTensor V,D;
    eigDecomp(Mr,V,D);

    CHECK((Mr*V - prime(V)*D).norm() < 1E-11);

    IQTensor M = Mr+Complex_i*Mi;
    eigDecomp(M,V,D);

    CHECK((M*V - prime(V)*D).norm() < 1E-11);
    }
}
