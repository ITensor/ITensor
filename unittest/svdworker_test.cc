#include "test.h"
#include "svdworker.h"
#include <boost/test/unit_test.hpp>

using namespace std;
using boost::format;

struct SVDWorkerDefaults
    {

    IQIndex S1,S2,L1,L2,Mid,Mid10;

    const Index
        s1u,s1d,s2u,s2d,
        l1u,l10,l1d,
        l2uu,l2u,l20,l2d,l2dd,
        mid,
        mid10u,mid10z,mid10d;

    IQTensor Phi0,L,R,Psi;
    IQTSparse V;
    ITensor phi0,l,r,psi;
    ITSparse v;

    SVDWorkerDefaults()
        :
        s1u(Index("Site1 Up",1,Site)),
        s1d(Index("Site1 Dn",1,Site)),
        s2u(Index("Site2 Up",1,Site)),
        s2d(Index("Site2 Dn",1,Site)),
        l1u(Index("Link1 Up",10,Link)),
        l10(Index("Link1 Z0",20,Link)),
        l1d(Index("Link1 Dn",10,Link)),
        l2uu(Index("Link2 UU",5,Link)),
        l2u(Index("Link2 Up",10,Link)),
        l20(Index("Link2 Z0",50,Link)),
        l2d(Index("Link2 Dn",10,Link)),
        l2dd(Index("Link2 DD",5,Link)),
        mid(Index("mid")),
        mid10u(Index("mid10u",2)),
        mid10z(Index("mid10z",5)),
        mid10d(Index("mid10d",3)) 
        {
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
                     l2uu,QN(+2),
                     l2u,QN(+1),
                     l20,QN( 0),
                     l2d,QN(-1),
                     l2dd,QN(-2),
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
            uudd.Randomize();
            Phi0 += uudd;

            ITensor zudz(l10,s1u,s2d,l20);
            zudz.Randomize();
            Phi0 += zudz;

            ITensor duud(l1d,s1u,s2u,l2d);
            duud.Randomize();
            Phi0 += duud;

            ITensor zuudd(l10,s1u,s2u,l2dd);
            zuudd.Randomize();
            Phi0 += zuudd;

            ITensor zdduu(l10,s1d,s2d,l2uu);
            zdduu.Randomize();
            Phi0 += zdduu;
            }
        Phi0 *= 1.0/Phi0.norm();
        checkDiv(Phi0);

        phi0 = Phi0;

        L = IQTensor(L1,S1,Mid10);
            {
            ITensor zud(l10,s1u,mid10d);
            zud.Randomize();
            L += zud;

            ITensor udz(l1u,s1d,mid10z);
            udz.Randomize();
            L += udz;

            ITensor duz(l1d,s1u,mid10z);
            duz.Randomize();
            L += duz;

            ITensor zdu(l10,s1d,mid10u);
            zdu.Randomize();
            L += zdu;
            }
        L *= 1./L.norm();
        l = L;

        R = IQTensor(Mid10,S2,L2);
            {
            ITensor zud(mid10z,s2u,l2d);
            zud.Randomize();
            R += zud;

            ITensor udz(mid10u,s2d,l20);
            udz.Randomize();
            R += udz;

            ITensor duz(mid10d,s2u,l20);
            duz.Randomize();
            R += duz;

            ITensor zdu(mid10z,s2d,l2u);
            zdu.Randomize();
            R += zdu;
            }
        R *= 1./R.norm();
        r = R;

        /*
        V = IQTSparse(Mid10);
            {
            Vector diag;

            diag.ReDimension(mid10u.m());
            for(int j = 1; j <= mid10u.m(); ++j)
                diag(j) = 0.5*sqrt(1.*j);
            V += ITSparse(mid10u,diag);

            diag.ReDimension(mid10z.m());
            for(int j = 1; j <= mid10z.m(); ++j)
                diag(j) = 2.121*sqrt(1.*j+7);
            V += ITSparse(mid10z,diag);

            diag.ReDimension(mid10d);
            for(int j = 1; j <= mid10d.m(); ++j)
                diag(j) = 4.323*sqrt(1.*j+10);
            V += ITSparse(mid10d,diag);
            }
        v = V;
        */

        Psi = L;
        //Psi /= V;
        Psi.conj(Mid10);
        Psi *= R;

        psi = Psi;
        }


    };

BOOST_FIXTURE_TEST_SUITE(SVDWorkerTest,SVDWorkerDefaults)

TEST(SiteSVD)
    {
    SVDWorker svd;

    //
    //ITensor version
    //

    ITensor a(L1,S1),b(S2,L2);

    //svd.showeigs(true);
    phi0 *= -1.2324;
    svd.denmatDecomp(phi0,a,b,Fromleft);

    //Print(((a*b)-phi0).norm());
    CHECK(((a*b)-phi0).norm() < 1E-12 );

    CHECK(svd.truncerr() < 1E-12);

    //
    //IQTensor version
    //

    IQTensor A(L1,S1),B(S2,L2);

    Phi0 *= -10.23432;
    svd.denmatDecomp(Phi0,A,B,Fromleft);

    //Print(((A*B)-Phi0).norm());
    CHECK(((A*B)-Phi0).norm() < 1E-12 );

    CHECK(svd.truncerr() < 1E-12);
    }

TEST(BondSVD)
    {
    SVDWorker svd;

    //
    //ITensor version
    //

    ITensor l(L1,S1,Mid),r(Mid,S2,L2);
    ITSparse v(Mid);

    phi0 *= -0.235;

    svd.csvd(phi0,l,v,r);

    //Print(((l*v*r)-phi0).norm());
    CHECK(((l*v*r)-phi0).norm() < 1E-12 );


    CHECK(svd.truncerr() < 1E-12);

    //
    //IQTensor version
    //

    IQTensor L(L1,S1,Mid),R(Mid,S2,L2);
    IQTSparse V(Mid);

    svd.csvd(Phi0,L,V,R);
    
    IQTensor nPhi = L*V*R;

    CHECK((nPhi-Phi0).norm() < 1E-12 );

    CHECK(svd.truncerr() < 1E-12);
    
    }

TEST(CSVDNorm)
    {
    SVDWorker svd;

    //
    //ITensor version
    //

    Index rr("rr",1,Link),s1("s1",3,Site),s2("s2",3,Site);
    phi0 = ITensor(s1,s2,rr);

    ITensor l(s1,Mid),r(Mid,s2,rr);
    ITSparse v(Mid);

    phi0(s1(3),s2(1),rr(1)) = -0.172148;
    phi0(s1(2),s2(2),rr(1)) = 0.427132;
    phi0(s1(1),s2(3),rr(1)) = -0.88765;
    phi0.scaleTo(-2.127145);

    svd.svd(phi0,l,v,r);
    Index nmr = index_in_common(r,v,Link);
    //ITensor rhoRsvd1 = r * conj(primeind(r,rr));
    //ITensor rhoRsvd2 = r * conj(primeind(r,nmr));
    //PrintDat(r);
    //PrintDat(rhoRsvd1);
    //PrintDat(rhoRsvd2);


    svd.csvd(phi0,l,v,r);

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

TEST(AbsoluteCutoff)
    {
    SVDWorker svd;

    svd.absoluteCutoff(true);

    //
    //ITensor version
    //

    ITensor a(L1,S1,Mid),b(Mid,S2,L2);
    ITSparse c(Mid);

    Real cutoff = 1E-3;
    svd.cutoff(cutoff);
    svd.csvd(phi0,a,c,b);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-5;
    svd.cutoff(cutoff);
    svd.csvd(phi0,a,c,b);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-7;
    svd.cutoff(cutoff);
    svd.csvd(phi0,a,c,b);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    //
    //IQTensor version
    //

    IQTensor A(L1,S1,Mid),B(Mid,S2,L2);
    IQTSparse C(Mid);

    cutoff = 1E-3;
    svd.cutoff(cutoff);
    svd.csvd(Phi0,A,C,B);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-5;
    svd.cutoff(cutoff);
    svd.csvd(Phi0,A,C,B);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-7;
    svd.cutoff(cutoff);
    svd.csvd(Phi0,A,C,B);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);
    }

/*
TEST(UseOrigM)
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
TEST(SvdRank2)
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
    uu.Randomize(); vv.Randomize();
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

    ITSparse D;
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
    IQTSparse DD;

    svd.svdRank2(TT,UU,DD,VV);

    }
    */


BOOST_AUTO_TEST_SUITE_END()
