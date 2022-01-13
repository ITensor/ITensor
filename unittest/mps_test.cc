#include "test.h"
#include "itensor/mps/mps.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/electron.h"
#include "itensor/mps/sites/fermion.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"
#include "itensor/util/iterate.h"
#include "mps_mpo_test_helper.h"
#include <iomanip>

using namespace itensor;
using std::vector;

//
//  os << Matrix1 is currently not supported in mat.h so we make a
//  temporary function to dump a Matrix1 object for viewing the output of the expect function.
//
std::ostream& operator<< (std::ostream& s, const Matrix1& M)
{
    for(auto r : range1(nrows(M)))
        {
        s << "|";
        for(auto c : range1(ncols(M)))
            {
            s << formatVal(M(r,c));
            s << (c == ncols(M) ? "|" : " ");
            }
        if(r < nrows(M)) s << "\n";
        }
        return s;
}


TEST_CASE("MPSTest")
{

static const int N = 10;

// SiteSet with no QNs
auto shsites = SpinHalf(N,{"ConserveQNs=",false});
auto shFerro = InitState(shsites,"Up");
auto shNeel = InitState(shsites);
for(auto j : range1(N))
    shNeel.set(j,j%2==1 ? "Up" : "Dn");


auto shsitesQNs = SpinHalf(N,{"ConserveQNs=",true});
auto shFerroQNs = InitState(shsitesQNs,"Up");
auto shNeelQNs = InitState(shsitesQNs);
for(auto j : range1(N))
    shNeelQNs.set(j,j%2==1 ? "Up" : "Dn");

SECTION("QNCheck")
    {
    auto psiNeelQNs = MPS(shNeelQNs);

    CHECK(checkQNs(psiNeelQNs));
    CHECK_EQUAL(totalQN(psiNeelQNs),QN({"Sz",0}));

    auto psiFerroQNs = MPS(shFerroQNs);

    CHECK(checkQNs(psiFerroQNs));
    CHECK_EQUAL(totalQN(psiFerroQNs),QN({"Sz",10}));
    }

SECTION("hasQNs")
    {
    CHECK(hasQNs(shFerroQNs));
    CHECK(not hasQNs(shFerro));
    }

SECTION("Constructors (m==1)")
    {
    auto psi = MPS(shsitesQNs);
    CHECK(1==dim(linkIndex(psi,2)));
    CHECK(checkTags(psi));
    }

SECTION("Constructors (dim>1)")
    {
    auto d = 4;
    auto psi = MPS(shsites,d);
    auto l2 = commonIndex(psi(2),psi(3));
    CHECK(d==dim(l2));

    CHECK(checkTags(psi));

    for(int n = 1; n <= N; ++n)
      psi.ref(n).randomize({"Complex=",true});

    psi.position(1);

    CHECK(checkTags(psi));
    CHECK(isOrtho(psi));

    psi.position(N);

    CHECK(checkTags(psi));
    CHECK(isOrtho(psi));

    psi.normalize();

    CHECK_CLOSE(norm(psi),1.);

    }

SECTION("Random constructors (dim==1)")
    {
    auto psi = randomMPS(shsites);
    psi.ref(1) *= Complex_i;

    CHECK(checkTags(psi));

    auto psi1 = psi;
    psi1.position(1);

    auto norm2 = innerC(psi1,psi1);

    CHECK_CLOSE(norm(psi1),std::sqrt(real(norm2)));
    CHECK_CLOSE(imag(norm2),0.);
    CHECK(checkTags(psi1));
    CHECK_CLOSE(diff(psi,psi1),0.);

    psi1.position(N);

    norm2 = innerC(psi1,psi1);

    CHECK_CLOSE(norm(psi1),std::sqrt(real(norm2)));
    CHECK_CLOSE(imag(norm2),0.);
    CHECK(checkTags(psi1));
    CHECK_CLOSE(diff(psi,psi1),0.);
    }

SECTION("Check .position() with custom tags")
    {
    auto psi = randomMPS(shNeelQNs);
    psi.ref(1) *= Complex_i;

    CHECK(checkTags(psi,"Site","Link"));

    psi.replaceTags("Site","MySite");
    psi.replaceTags("Link","MyLink");

    CHECK(checkTags(psi,"MySite","MyLink"));

    auto opsi = psi;

    opsi.position(1);

    auto norm_opsi = std::sqrt(real(innerC(opsi,opsi)));

    CHECK(checkOrtho(opsi));
    CHECK_CLOSE(norm(opsi),norm_opsi);
    CHECK_CLOSE(diff(opsi,psi),0.);
    CHECK(checkTags(opsi,"MySite","MyLink"));

    opsi.position(N);

    CHECK_CLOSE(diff(opsi,psi),0.);
    CHECK(checkTags(opsi,"MySite","MyLink"));
    }

SECTION("Check .orthogonalize() with custom tags")
    {
    auto psi = randomMPS(shNeelQNs);
    psi.ref(1) *= Complex_i;

    CHECK(checkTags(psi,"Site","Link"));

    psi.replaceTags("Site","MySite");
    psi.replaceTags("Link","MyLink");

    CHECK(checkTags(psi,"MySite","MyLink"));

    auto opsi = psi;

    opsi.orthogonalize();

    CHECK(checkOrtho(opsi));
    CHECK_CLOSE(norm(opsi),std::sqrt(real(innerC(opsi,opsi))));
    CHECK_CLOSE(diff(opsi,psi),0.);
    CHECK(checkTags(opsi,"MySite","MyLink"));
    }

SECTION("Random constructors, QN conserved (dim==1)")
    {
    auto psi = randomMPS(shFerroQNs);
    psi.ref(1) *= Complex_i;
    psi.orthogonalize();
    CHECK(norm(psi)>0);
    CHECK(checkTags(psi));
    }

SECTION("MPSAddition 1")
    {
    auto sites = Fermion(10);

    auto i1 = InitState(sites,"Emp");
    auto i2 = i1;

    i1.set(1,"Occ");
    i2.set(2,"Occ");

    //"Valence bond" between sites 1 and 2
    auto psi = ISqrt2*sum(MPS(i1),MPS(i2));

    CHECK(checkTags(psi));
    CHECK_CLOSE(norm(psi),1);
    CHECK_EQUAL(totalQN(psi),QN({"Nf",1}));
    }

SECTION("MPSAddition Custom Tags")
    {
    auto sites = Fermion(10,{"ConserveQNs=",true});

    auto i1 = InitState(sites,"Emp");
    auto i2 = i1;

    i1.set(1,"Occ");
    i2.set(2,"Occ");

    auto psi1 = MPS(i1);
    psi1.ref(1) *= Complex_i;
    auto psi2 = MPS(i2);
    psi2.ref(1) *= Complex_i;

    psi1.replaceTags("Site","MySite").replaceTags("Link","MyLink1");
    psi2.replaceTags("Site","MySite").replaceTags("Link","MyLink2");

    //"Valence bond" between sites 1 and 2
    auto psi12 = ISqrt2*sum(psi1,psi2);

    CHECK(checkTags(psi12,"MySite","MyLink1"));
    CHECK_CLOSE(norm(psi12),1);
    CHECK_EQUAL(totalQN(psi12),QN({"Nf",1}));

    auto psi21 = ISqrt2*sum(psi2,psi1);
    CHECK(checkTags(psi21,"MySite","MyLink2"));
    CHECK_CLOSE(norm(psi21),1);
    CHECK_EQUAL(totalQN(psi21),QN({"Nf",1}));
    }

SECTION("MPSAddition 2")
    {
    auto sites = SiteSet(10,2);
    auto psi1 = randomMPS(sites);
    psi1.ref(1) *= Complex_i;
    auto psi2 = randomMPS(sites);
    psi2.ref(1) *= Complex_i;
    auto psi = sum(psi1,psi2);

    auto exact_inner = innerC(psi1,psi2);
    exact_inner += innerC(psi2,psi1);
    exact_inner += innerC(psi1,psi1);
    exact_inner += innerC(psi2,psi2);

    CHECK(checkTags(psi));
    CHECK_CLOSE(exact_inner,innerC(psi,psi));
    CHECK_EQUAL(order(psi(1)),2);
    CHECK_EQUAL(order(psi(2)),3);
    CHECK_EQUAL(order(psi(5)),3);
    CHECK_EQUAL(order(psi(9)),3);
    CHECK_EQUAL(order(psi(10)),2);
    }

SECTION("PositionTest")
    {
    auto sites = Fermion(10);

    auto init = InitState(sites,"Emp");
    init.set(2,"Occ");
    init.set(4,"Occ");
    init.set(6,"Occ");

    auto psi = MPS(init);
    psi.ref(1) *= Complex_i;

    psi.position(1,{"Cutoff=",1E-8});

    CHECK(checkOrtho(psi));
    CHECK_EQUAL(findCenter(psi),1);

    psi.position(4,{"Cutoff=",1E-8});

    CHECK_EQUAL(findCenter(psi),4);
    }

SECTION("Orthogonalize")
    {
    auto d = 20;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    //Make a random MPS of bond dim. d
    auto psi = MPS(sites,d);
    for(auto n : range1(N))
        psi.ref(n).randomize({"Complex=",true});

    CHECK(checkTags(psi));

    //Normalize psi
    auto n2 = real(innerC(psi,psi));
    psi.ref(1) /= sqrt(n2);

    auto opsi = psi;
    opsi.orthogonalize({"Cutoff",1E-16});

    CHECK(checkOrtho(opsi));
    CHECK_CLOSE(innerC(opsi,psi),1.0);
    CHECK(checkTags(opsi));

    psi.orthogonalize({"MaxDim=",10,"Cutoff=",1E-16});

    CHECK(checkOrtho(psi));
    CHECK(maxLinkDim(psi)==10);
    CHECK(checkTags(psi));
    }

SECTION("Overlap - 1 site")
    {
    auto psi = MPS(1);
    auto s = Index(2,"s");
    psi.ref(1) = randomITensorC(s);

    CHECK_CLOSE(innerC(psi,psi),eltC(dag(psi(1))*psi(1)));
    }

SECTION("siteInds")
    {
    auto sites = SpinHalf(N);
    auto initstate = InitState(sites,"Up");
    auto psi = randomMPS(initstate);

    CHECK(checkTags(psi));

    auto s = siteInds(psi);

    for( auto n : range1(N) )
      {
      CHECK( siteIndex(psi,n)==s(n) );
      CHECK( sites(n)==s(n) );
      }
    }

SECTION("replaceSiteInds")
    {
    auto sites1 = SpinHalf(N);
    auto initstate1 = InitState(sites1,"Up");
    auto psi1 = randomMPS(initstate1);
    psi1.position(1);
    psi1 /= norm(psi1);

    CHECK(checkTags(psi1));

    auto sites2 = SpinHalf(N);
    auto initstate2 = InitState(sites2,"Up");
    auto psi2 = randomMPS(initstate2);
    psi2.position(1);
    psi2 /= norm(psi2);
    psi2.replaceTags("Site","MySite");

    CHECK(checkTags(psi2,"MySite","Link"));

    auto s1 = siteInds(psi1);
    auto s2 = siteInds(psi2);

    auto psi1_new = replaceSiteInds(psi1,s2);

    CHECK(checkTags(psi1_new,"MySite","Link"));
    for( auto n : range1(N) )
      CHECK( siteIndex(psi1_new,n)==siteIndex(psi2,n) );
    }

SECTION("prime")
    {
    auto s = SpinHalf(N,{"ConserveQNs=",false});
    auto psi = randomMPS(s);

    auto psi1 = prime(psi);

    CHECK( prime(siteIndex(psi,3)) == siteIndex(psi1,3) );
    CHECK( prime(linkIndex(psi,3)) == linkIndex(psi1,3) );

    auto psi2 = prime(psi,"Site");

    CHECK( prime(siteIndex(psi,3)) == siteIndex(psi2,3) );
    CHECK( linkIndex(psi,3) == linkIndex(psi2,3) );

    auto psi3 = mapPrime(psi2,1,2);

    CHECK( prime(siteIndex(psi,3),2) == siteIndex(psi3,3) );
    CHECK( linkIndex(psi,3) == linkIndex(psi3,3) );

    auto psi4 = swapPrime(psi2,0,1);
    
    CHECK( siteIndex(psi,3) == siteIndex(psi4,3) );
    CHECK( prime(linkIndex(psi,3)) == linkIndex(psi4,3) );

    }

//-----------------------------------------------------------------------------------
//
//  Begin testing for the expect(psi,ops) function.
//
SECTION("expect S1/2 No QNs ferro")
{
    SiteSet   sites = SpinHalf(N, {"ConserveQNs=",false});
    MPS       psi   = randomMPS(InitState(sites,"Up"));

    // Test with the site range option
    int site_start=2,site_end=4;
    Matrix1 ex=expect(psi, {site_start,site_end},"Sz","ISy","Sx","S+","S-","S2","Sz*Sz","ISy*ISy","Sx*Sx","projUp","projDn") ;

//    std::cout << "expect table" << std::endl;
//    std::cout << "   Sz        ISy       Sx        S+        S-        S2        Sz*Sz     Sy*Sy     Sx*Sx     projUp    projDn" << std::endl;
//    std::cout << ex << std::endl;

    for (int j=1; j<=site_end-site_start+1; j++) //site loop
    {
        double Sz=ex(j,1),ISy=ex(j,2),Sx=ex(j,3),Sp=ex(j,4),Sm=ex(j,5),S2=ex(j,6);
        CHECK_CLOSE(sqrt(Sz*Sz+ISy*ISy+Sx*Sx),0.5);
        CHECK_CLOSE(sqrt(Sz*Sz+0.5*(Sm*Sp+Sp*Sm)),0.5);
        CHECK_CLOSE(S2,0.75);
        CHECK_CLOSE(Sx,0.5*(Sp+Sm));
        CHECK_CLOSE(ISy,0.5*(Sp-Sm));
    }
}

SECTION("expect S1/2 With QNs ferro")
{
    SiteSet   sites = SpinHalf(N);
    MPS       psi   = randomMPS(InitState(sites,"Up"));

    // Only ops that commute with Sz are allowed here.
    Matrix1 ex=expect(psi,"Sz","S+","S-","S2","Sz*Sz","projUp","projDn") ;

//    std::cout << "expect table" << std::endl;
//    std::cout << "  Sz        S+        S-        S2        Sz*Sz     projUp    projDn" << std::endl;
//    std::cout << ex << std::endl;

    for (int j=1; j<=psi.length(); j++) //site loop
    {
        double Sz=ex(j,1),Sp=ex(j,2),Sm=ex(j,3),S2=ex(j,4),SzSz=ex(j,5),projUp=ex(j,6),projDn=ex(j,7);
        CHECK_CLOSE(Sz,0.5);
        CHECK_CLOSE(SzSz,0.25);
        CHECK_CLOSE(Sp,0.0);
        CHECK_CLOSE(Sm,0.0);
        CHECK_CLOSE(projUp,1.0);
        CHECK_CLOSE(projDn,0.0);
        CHECK_CLOSE(S2,0.75);
    }
}

SECTION("expect Electron With QNs ")
{
    SiteSet   sites = Electron(N);
    MPS       psi   = randomMPS(InitState(sites,"Up"));

    Matrix1 ex=expect(psi,"Sz","S+","S-","S2","Nup","Ndn","Nupdn","Ntot") ;

//    std::cout << "expect table" << std::endl;
//    std::cout << "  Sz        S+        S-        S2        Nup       Ndn       NupDn     Ntot" << std::endl;
//    std::cout << ex << std::endl;

    for (int j=1; j<=psi.length(); j++) //site loop
    {
        double Sz=ex(j,1),Sp=ex(j,2),Sm=ex(j,3),S2=ex(j,4),Nup=ex(j,5),Ndn=ex(j,6),Nupdn=ex(j,7),Ntot=ex(j,8);
        CHECK_CLOSE(sqrt(Sz*Sz+0.5*(Sm*Sp+Sp*Sm)),0.5);
        CHECK_CLOSE(S2,0.75);
        CHECK_CLOSE(Nup,1.0);
        CHECK_CLOSE(Ndn,0.0);
        CHECK_CLOSE(Nupdn,0.0);
        CHECK_CLOSE(Ntot,1.0);
    }
}


SECTION("expect Fermion No QNs ")
{
    SiteSet   sites = Fermion(N, {"ConserveQNs=",false});
    MPS       psi   = randomMPS(sites);

    Matrix1 ex=expect(psi,"N","Cdag*C","Adag*A","F","projEmp","projOcc") ;

//    std::cout << "expect table" << std::endl;
//    std::cout << "    N        Cdag*C      Adag*A     F     projEmp   projOcc" << std::endl;
//    std::cout << ex << std::endl;

    for (int j=1; j<=psi.length(); j++) //site loop
    {
        double NN=ex(j,1),CdagC=ex(j,2),AdagA=ex(j,3),F=ex(j,4),projEmp=ex(j,5),projOcc=ex(j,6);
        CHECK_CLOSE(F,1-2*NN);
        CHECK_CLOSE(CdagC,NN);
        CHECK_CLOSE(AdagA,NN);
        CHECK_CLOSE(projOcc,NN);
        CHECK_CLOSE(projOcc+projEmp,1.0);
    }
}

SECTION("expect Fermion With QNs ")
{
    SiteSet   sites = Fermion(N, {"ConserveQNs=",true});
    MPS       psi   = randomMPS(InitState(sites,"1"));

    Matrix1 ex=expect(psi,"N","Cdag*C","Adag*A","F","projEmp","projOcc") ;

//    std::cout << "expect table" << std::endl;
//    std::cout << "    N       Cdag*C      Adag*A     F      projEmp   projOcc" << std::endl;
//    std::cout << ex << std::endl;

    for (int j=1; j<=psi.length(); j++) //site loop
    {
        double NN=ex(j,1),CdagC=ex(j,2),AdagA=ex(j,3),F=ex(j,4),projEmp=ex(j,5),projOcc=ex(j,6);
        CHECK_CLOSE(NN,1.0);
        CHECK_CLOSE(F,1.0-2*NN);
        CHECK_CLOSE(CdagC,NN);
        CHECK_CLOSE(AdagA,NN);
        CHECK_CLOSE(projOcc,NN);
        CHECK_CLOSE(projOcc+projEmp,1.0);
    }
}

}
