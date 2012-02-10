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

    IQTensor Phi0,L,V,R,Psi;
    ITensor phi0,l,v,r,psi;

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
        Phi0.checkDiv();

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

        V = IQTensor(Mid10);
            {
            ITensor u(mid10u);
            for(int j = 1; j <= mid10u.m(); ++j)
                u(mid10u(j)) = 0.5*sqrt(1.*j);
            V += u;

            ITensor z(mid10z);
            for(int j = 1; j <= mid10z.m(); ++j)
                z(mid10z(j)) = 2.121*sqrt(1.*j+7);
            V += z;

            ITensor d(mid10d);
            for(int j = 1; j <= mid10d.m(); ++j)
                d(mid10d(j)) = 4.323*sqrt(1.*j+10);
            V += d;
            }
        v = V;

        Psi = L/V;
        Psi.conj(V.index(1));
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
    svd(phi0,a,b,Fromleft);

    Print(((a*b)-phi0).norm());
    CHECK(((a*b)-phi0).norm() < 1E-12 );

    CHECK(svd.truncerr() < 1E-12);

    //
    //IQTensor version
    //

    IQTensor A(L1,S1),B(S2,L2);

    svd(Phi0,A,B,Fromleft);

    CHECK(((A*B)-Phi0).norm() < 1E-12 );

    CHECK(svd.truncerr() < 1E-12);
    }

TEST(BondSVD)
    {
    SVDWorker svd;

    //
    //ITensor version
    //

    ITensor l(L1,S1,Mid),r(Mid,S2,L2),v(Mid);

    svd(phi0,l,v,r);

    CHECK((((l/v)*r)-phi0).norm() < 1E-12 );

    CHECK(((l*(v/r))-phi0).norm() < 1E-12 );

    CHECK(svd.truncerr() < 1E-12);

    //
    //IQTensor version
    //

    IQTensor L(L1,S1,Mid),R(Mid,S2,L2),V(Mid);

    svd(Phi0,L,V,R);
    
    IQTensor nPhi = L/V;
    nPhi.conj(V.index(1));
    nPhi *= R;

    CHECK((nPhi-Phi0).norm() < 1E-12 );

    CHECK(svd.truncerr() < 1E-12);
    }

TEST(AbsoluteCutoff)
    {
    SVDWorker svd;

    svd.absoluteCutoff(true);

    //
    //ITensor version
    //

    ITensor a(L1,S1,Mid),c(Mid),b(Mid,S2,L2);

    Real cutoff = 1E-3;
    svd.cutoff(cutoff);
    svd(phi0,a,c,b);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-5;
    svd.cutoff(cutoff);
    svd(phi0,a,c,b);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-7;
    svd.cutoff(cutoff);
    svd(phi0,a,c,b);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    //
    //IQTensor version
    //

    IQTensor A(L1,S1,Mid),C(Mid),B(Mid,S2,L2);

    cutoff = 1E-3;
    svd.cutoff(cutoff);
    svd(Phi0,A,C,B);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-5;
    svd.cutoff(cutoff);
    svd(Phi0,A,C,B);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);

    cutoff = 1E-7;
    svd.cutoff(cutoff);
    svd(Phi0,A,C,B);
    CHECK(svd.eigsKept()(svd.numEigsKept()) > cutoff);
    }

TEST(Minm)
    {
    SVDWorker svd;

    svd.cutoff(1E-3);

    //
    //ITensor version
    //

    ITensor a(L1,S1,Mid),c(Mid),b(Mid,S2,L2);

    svd.minm(50);
    svd(phi0,a,c,b);

    CHECK_EQUAL(50,c.index(1).m());
        
    svd.minm(75);
    svd(phi0,a,c,b);

    CHECK_EQUAL(75,c.index(1).m());

    //Check that if the requested minm is larger
    //than the number of possibly non-zero 
    //denmat eigs, we just keep all of them
    const int max_possible = min(L1.m()*S1.m(),S2.m()*L2.m());

    svd.minm(20*max_possible);
    svd(phi0,a,c,b);
    Print(c.index(1).m());

    //Commented out until I fix the way QNs work
    //in canonical MPS
    //CHECK_EQUAL(max_possible,c.index(1).m());

    //
    //IQTensor version
    //

    IQTensor A(L1,S1,Mid),C(Mid),B(Mid,S2,L2);

    svd.minm(50);
    svd(Phi0,A,C,B);

    CHECK_EQUAL(50,C.index(1).m());

    const int iq_max_possible = 70;
        
    svd.minm(20*iq_max_possible);
    svd(Phi0,A,C,B);

    //Commented out until I fix the way QNs work
    //in canonical MPS
    //CHECK_EQUAL(iq_max_possible,C.index(1).m());

    IQTensor nPhi = A/C;
    nPhi.conj(C.index(1));
    nPhi *= B;
    Print((nPhi-Phi0).norm());
    CHECK((nPhi-Phi0).norm() < 1E-12);

    }

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

    svd(Psi,L,V,R);

    CHECK_EQUAL(Mid10.m(),V.index(1).m());

    IQTensor nPsi = (L/V);
    nPsi.conj(V.index(1));
    nPsi *= R;

    CHECK((nPsi-Psi).norm() < 1E-12);
    CHECK(svd.truncerr() < 1E-12);

    //
    //ITensor version
    //

    svd(psi,l,v,r);

    CHECK_EQUAL(Mid10.m(),v.index(1).m());
    
    CHECK((((l/v)*r)-psi).norm() < 1E-12 );
    CHECK(svd.truncerr() < 1E-12);
    }


BOOST_AUTO_TEST_SUITE_END()
