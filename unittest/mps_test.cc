#include "test.h"
#include "itensor/mps/mps.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/electron.h"
#include "itensor/mps/sites/fermion.h"
#include "itensor/mps/autompo.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"
#include "itensor/util/iterate.h"
#include "mps_mpo_test_helper.h"
#include <iomanip>

using namespace itensor;
using std::vector;

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

SECTION("Random with prescribed bond dimensions")
    {
    auto psi = randomMPS(shsites,1);
    for(auto n : range1(N-1))
      {
      CHECK_EQUAL(1,dim(linkIndex(psi,n)));
      }
    psi = randomMPS(shsites,2);
    for(auto n : range1(N-1))
      {
      CHECK_EQUAL(2,dim(linkIndex(psi,n)));
      }
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
} //TEST_CASE("MPSTest")

//-----------------------------------------------------------------------------------
//
//  Begin testing for the expect(psi,ops) function.
//
//  un-comment this define if you want to see output of tables from the expect tests
//#define EXPECT_VERBOSE

//
//  Helper functions
//    os << VecT and os << VecVecT for viewing the output of the expect 
//    and correlationMatrix functions.
//

template <typename T>
std::ostream& operator<< (std::ostream& s, const vector<T>& v)
{
    s << "|";
    for(auto& c : v)
        {
        if (std::fabs(std::imag(c))<1e-15)
            s << formatVal(std::real(c));
        else
            s << formatVal(c);
            
        s << (&c == &v.back() ? "|" : " ");
        }
    s << "\n";
    return s;
}

template <typename T>
std::ostream& operator<< (std::ostream& s, const vector<vector<T>>& m)
{
    for(auto& r : m)
        s << r;
    return s;
}

TEST_CASE("expect function")
{
    int N=10;

    SECTION("expect Real Psi, Real Ops, S1/2 No QNs, range")
    {
        SiteSet   sites = SpinHalf(N, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up")); 

        // Test with the site range option
        auto ex=expect(psi,sites,{"Sz","ISy","Sx","S+","S-","S2","Sz*Sz","ISy*ISy","Sx*Sx","projUp","projDn"},range1(2,5)) ;
#ifdef EXPECT_VERBOSE
        std::cout << "expect table" << std::endl;
        std::cout << "   Sz        ISy       Sx        S+        S-        S2        Sz*Sz    ISy*ISy    Sx*Sx     projUp    projDn" << std::endl;
        std::cout << ex << std::endl;
#endif
        for (auto& j:ex) //site loop
        {
            auto i=j.begin(); //FYI type of i should be: VecVecR::value_type::const_iterator
            Real Sz=*i++,ISy=*i++,Sx=*i++,Sp=*i++,Sm=*i++,S2=*i++;
            CHECK_CLOSE(sqrt(Sz*Sz+ISy*ISy+Sx*Sx),0.5);
            CHECK_CLOSE(sqrt(Sz*Sz+0.5*(Sm*Sp+Sp*Sm)),0.5);
            CHECK_CLOSE(S2,0.75);
            CHECK_CLOSE(Sx,0.5*(Sp+Sm));
            CHECK_CLOSE(ISy,0.5*(Sp-Sm));
        }
        
        // Test single operator version
        auto ex1=expect(psi,sites,"Sz",range1(2,5)) ;
        for (VecR::size_type i=0;i<ex1.size();i++)
            CHECK(ex1[i]==ex[i][0]);
        
    }

    SECTION("expect Real Psi, Complex Ops, S1/2 No QNs, range")
    {
        SiteSet   sites = SpinHalf(N, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up")); 

        // Test with the site range option
        auto ex=expectC(psi,sites,{"Sz","Sy","Sx","S+","S-","S2","Sz*Sz","Sy*Sy","Sx*Sx","projUp","projDn"},range1(2,5)) ;
        // ex should be of type: VecVecC
#ifdef EXPECT_VERBOSE
        std::cout << "expect table" << std::endl;
        std::cout << "   Sz        Sy       Sx        S+        S-        S2        Sz*Sz     Sy*Sy     Sx*Sx     projUp    projDn" << std::endl;
        std::cout << ex << std::endl;
#endif
        for (auto& j:ex) //site loop
        {
            auto i=j.begin();
            Complex Sz=*i++,Sy=*i++,Sx=*i++,Sp=*i++,Sm=*i++,S2=*i++;
            CHECK_CLOSE(sqrt(Sz*Sz+Sy*Sy+Sx*Sx),0.5);
            CHECK_CLOSE(sqrt(Sz*Sz+0.5*(Sm*Sp+Sp*Sm)),0.5);
            CHECK_CLOSE(S2,0.75);
            CHECK_CLOSE(Sx,0.5*(Sp+Sm));
            CHECK_CLOSE(Sy,0.5*(Sp-Sm));
        }
        // Test single operator version
        auto ex1=expectC(psi,sites,"Sz",range1(2,5)) ;
        for (VecR::size_type i=0;i<ex1.size();i++)
            CHECK(ex1[i]==ex[i][0]);
    }

    SECTION("expect Complex Psi, Complex Ops, S1/2 No QNs, range")
    {
        SiteSet   sites = SpinHalf(N, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up"),{"Complex=",true}); 

        // Test with the site range option
        auto ex=expectC(psi,sites,{"Sz","Sy","Sx","S+","S-","S2","Sz*Sz","Sy*Sy","Sx*Sx","projUp","projDn"},range1(2,5)) ;
#ifdef EXPECT_VERBOSE
        std::cout << "expect table" << std::endl;
        std::cout << "      Sz        Sy       Sx             S+                 S-        S2        Sz*Sz     Sy*Sy     Sx*Sx     projUp    projDn" << std::endl;
        std::cout << ex << std::endl;
#endif
        for (auto& j:ex) //site loop
        {
            auto i=j.begin();
            Complex Sz=*i++,Sy=*i++,Sx=*i++,Sp=*i++,Sm=*i++,S2=*i++;
            CHECK_CLOSE(sqrt(Sz*Sz+Sy*Sy+Sx*Sx),0.5);
            CHECK_CLOSE(sqrt(Sz*Sz+0.5*(Sm*Sp+Sp*Sm)),0.5);
            CHECK_CLOSE(S2,0.75);
            CHECK_CLOSE(Sx,0.5*(Sp+Sm));
            CHECK_CLOSE(Sy,Complex(0,-0.5)*(Sp-Sm));
        }
        // Test single operator version
        auto ex1=expectC(psi,sites,"Sz",range1(2,5)) ;
        for (VecR::size_type i=0;i<ex1.size();i++)
            CHECK(ex1[i]==ex[i][0]);
    }


    SECTION("expect Real S1/2 With QNs ferro, site list")
    {
        SiteSet   sites = SpinHalf(N);
        MPS       psi   = randomMPS(InitState(sites,"Up"));

        // Only ops that commute with Sz are allowed here.
        auto ex=expect(psi,sites,{"Sz","S+","S-","S2","Sz*Sz","projUp","projDn"},{1,3,5,7,9}) ;
#ifdef EXPECT_VERBOSE
        std::cout << "expect table" << std::endl;
        std::cout << "  Sz        S+        S-        S2        Sz*Sz     projUp    projDn" << std::endl;
        std::cout << ex << std::endl;
#endif
        for (auto& j:ex) //site loop
        {
            auto i=j.begin();
            Real Sz=*i++,Sp=*i++,Sm=*i++,S2=*i++,SzSz=*i++,projUp=*i++,projDn=*i++;
            CHECK_CLOSE(Sz,0.5);
            CHECK_CLOSE(SzSz,0.25);
            CHECK_CLOSE(Sp,0.0);
            CHECK_CLOSE(Sm,0.0);
            CHECK_CLOSE(projUp,1.0);
            CHECK_CLOSE(projDn,0.0);
            CHECK_CLOSE(S2,0.75);
        }
        // Test single operator version
        auto ex1=expect(psi,sites,"Sz",{1,3,5,7,9}) ;
        for (VecR::size_type i=0;i<ex1.size();i++)
            CHECK(ex1[i]==ex[i][0]);
    }


    SECTION("expect Electron With QNs, no range, no site list ")
    {
        SiteSet   sites = Electron(N);
        MPS       psi   = randomMPS(InitState(sites,"Up"));

        auto ex=expect(psi,sites,{"Sz","S+","S-","S2","Nup","Ndn","Nupdn","Ntot"}) ;
#ifdef EXPECT_VERBOSE
        std::cout << "expect table" << std::endl;
        std::cout << "  Sz        S+        S-        S2        Nup       Ndn       NupDn     Ntot" << std::endl;
        std::cout << ex << std::endl;
#endif
        for (auto& j:ex) //site loop
        {
            auto i=j.begin();
            Real Sz=*i++,Sp=*i++,Sm=*i++,S2=*i++,Nup=*i++,Ndn=*i++,Nupdn=*i++,Ntot=*i++;
            CHECK_CLOSE(sqrt(Sz*Sz+0.5*(Sm*Sp+Sp*Sm)),0.5);
            CHECK_CLOSE(S2,0.75);
            CHECK_CLOSE(Nup,1.0);
            CHECK_CLOSE(Ndn,0.0);
            CHECK_CLOSE(Nupdn,0.0);
            CHECK_CLOSE(Ntot,1.0);
        }
        // Test single operator version
        auto ex1=expect(psi,sites,"Sz") ;
        for (VecR::size_type i=0;i<ex1.size();i++)
            CHECK(ex1[i]==ex[i][0]);
    }


    SECTION("expect Fermion No QNs ")
    {
        SiteSet   sites = Fermion(N, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(sites);

        auto ex=expect(psi,sites,{"N","Cdag*C","Adag*A","F","projEmp","projOcc"}) ;
#ifdef EXPECT_VERBOSE
        std::cout << "expect table" << std::endl;
        std::cout << "    N        Cdag*C      Adag*A     F     projEmp   projOcc" << std::endl;
        std::cout << ex << std::endl;
#endif
        for (auto& j:ex) //site loop
        {
            auto i=j.begin();
            double NN=*i++,CdagC=*i++,AdagA=*i++,F=*i++,projEmp=*i++,projOcc=*i++;
            CHECK_CLOSE(F,1-2*NN);
            CHECK_CLOSE(CdagC,NN);
            CHECK_CLOSE(AdagA,NN);
            CHECK_CLOSE(projOcc,NN);
            CHECK_CLOSE(projOcc+projEmp,1.0);
        }
        // Test single operator version
        auto ex1=expect(psi,sites,"Cdag*C") ;
        for (VecR::size_type i=0;i<ex1.size();i++)
            CHECK(ex1[i]==ex[i][1]);
    }


    SECTION("expect Fermion With QNs ")
    {
        SiteSet   sites = Fermion(N, {"ConserveQNs=",true});
        MPS       psi   = randomMPS(InitState(sites,"1"));

        auto ex=expect(psi,sites,{"N","Cdag*C","Adag*A","F","projEmp","projOcc"}) ;
#ifdef EXPECT_VERBOSE
        std::cout << "expect table" << std::endl;
        std::cout << "    N       Cdag*C      Adag*A     F      projEmp   projOcc" << std::endl;
        std::cout << ex << std::endl;
#endif
        for (auto& j:ex) //site loop
        {
            auto i=j.begin();
            double NN=*i++,CdagC=*i++,AdagA=*i++,F=*i++,projEmp=*i++,projOcc=*i++;
            CHECK_CLOSE(NN,1.0);
            CHECK_CLOSE(F,1.0-2*NN);
            CHECK_CLOSE(CdagC,NN);
            CHECK_CLOSE(AdagA,NN);
            CHECK_CLOSE(projOcc,NN);
            CHECK_CLOSE(projOcc+projEmp,1.0);
        }
        // Test single operator version
        auto ex1=expect(psi,sites,"Cdag*C") ;
        for (VecR::size_type i=0;i<ex1.size();i++)
            CHECK(ex1[i]==ex[i][1]);
    }

} //TEST_CASE("expect function")

//-----------------------------------------------------------------------------------
//
//  Begin testing for the correlation_matrix(psi,sites,op1,op2,site_range) function.
//


//
//  Helper functions
//    Easy check: make sure the diagonal agree with expect.
//
template <typename T> 
void checkDiagonalWithExpect(const MPS& psi,
                   const SiteSet& sites,
                   const std::pair<string,string>& op,
                   const vector<vector<T>>& cm,
                   detail::RangeHelper<int> site_range=range1(0))
{
    vector<string> op12;
    op12.push_back(op.first+"*"+op.second);
    vector<vector<T>> ex=expectT<T>(psi,sites,op12,site_range);
    for (VecR::size_type i=0;i<cm.size();i++)
        CHECK_CLOSE(cm[i][i],ex[i][0]);    
}

//
//  Alternate version of the corr matrix calculated through autoMPO
//
template <typename T> 
vector<vector<T>>
AutoMPOCorrelationMatrix
( MPS& psi,
 const SiteSet& sites,
 const std::pair<string,string>& op,
 detail::RangeHelper<int> site_range=range1(0))
{
    fixRange(site_range,sites.length());
    psi.orthogonalize();
    psi.normalize();
    vector<vector<T>> cm;
    for (auto i:site_range)
    {
        vector<T> row;
        for(auto j:site_range)
        {
          auto a = AutoMPO(sites);
          a += op.first, i, op.second, j;
          T cAutoMPO=innerT<T>(psi, toMPO(a), psi);
          row.push_back(cAutoMPO);
        }
        cm.push_back(row);
    }
    return cm;   
}
//
//  Check two matrices are equal within tolerance of CHECK_NEAR
//
template <typename T> 
void checkMatrices
   (const vector<vector<T>>& cm,
    const vector<vector<T>>& cm_fromAutoMPO,
    const std::pair<string,string>& op)
{
    bool pass=true;
    CHECK_EQUAL(cm.size(),cm_fromAutoMPO.size());
    typename vector<vector<T>>::const_iterator row_auto=cm_fromAutoMPO.begin();
    for (auto row:cm)
    {
        CHECK_EQUAL(row.size(),row_auto->size());
        typename vector<T>::const_iterator c_auto=row_auto->begin();
        for(auto c:row)
        {
          CHECK_CLOSE(*c_auto,c);
          pass = pass && (std::norm(*c_auto-c) < 1E-10);
          c_auto++;
        }
        row_auto++;
    }
    if (!pass)
    {
        std::cout << "correlation matrix for <" << op.first << "_i*" << op.second << "_j>"  
            << std::endl << cm << std::endl;
        std::cout << "auto mpo    matrix for <" << op.first << "_i*" << op.second << "_j>"  
            << std::endl << cm_fromAutoMPO << std::endl;            
    }    
}



TEST_CASE("correlationMatrix function")
{
    int N=10,Nsmall=3; //Use small lattices since checks using autoMPO are CPU intensive. 

    SECTION("correlationMatrix Real Psi, Real Ops, S1/2 No QNs, range")
    {
        SiteSet   sites = SpinHalf(N, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up")); 
        
        // Test with the site range option
        auto r=range1(2,7); //Don't use range here!!
        auto Ns=last(r)-current(r); //Number of sites in the range
        assert(Ns==6);
        vector<std::pair<string,string>> ops({{"Sz","Sz"},{"S+","S-"},{"S-","S+"},{"ISy","ISy"}});

        for (auto op:ops)
        {
            auto cm=correlationMatrix(psi,sites,op.first,op.second,r) ;
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op,r);
            CHECK(cm.size()==Ns);
            CHECK(cm[0].size()==Ns);
            checkDiagonalWithExpect<Real>(psi,sites,op,cm,r);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
        
    }

    // This type of cross correlation makes an non-symmetric correlation matrix.
    SECTION("correlation_matrix Real Psi, Sz*Sx, S1/2 No QNs")
    {
        SiteSet   sites = SpinHalf(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up")); 
        
        vector<std::pair<string,string>> ops({{"Sx","Sz"}});

        for (auto op:ops)
        {
            string op12=op.first+"*"+op.second; //check accept "Sx*Sz" format
            auto cm        =       correlationMatrix      (psi,sites,op.first,op.second,{"isHermitian",false}) ;
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op);

            checkDiagonalWithExpect<Real>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op); 
        } //ops loop
        
    }


    SECTION("correlation_matrix Complex Psi, Complex Ops, S1/2 No QNs, range")
    {
        SiteSet   sites = SpinHalf(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up"),{"Complex=",true}); 
        
        vector<std::pair<string,string>> ops({{"Sz","Sz"},{"S+","S-"},{"S-","S+"},{"Sy","Sy"}});

        for (auto op:ops)
        {
            auto cm=               correlationMatrixC        (psi,sites,op.first,op.second) ;
            auto cm_autompo=AutoMPOCorrelationMatrix<Complex>(psi,sites,op);

            checkDiagonalWithExpect<Complex>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
        
    }

    // This type of cross correlation makes a non-symmetric correlation matrix.
    SECTION("correlation_matrix Complex Psi, Sz*Sx, S1/2 No QNs, isHermitian=false ")
    {
        SiteSet   sites = SpinHalf(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up"),{"Complex=",true}); 
        
        vector<std::pair<string,string>> ops({{"Sx","Sz"}});

        for (auto op:ops)
        {
            auto cm=correlationMatrixC(psi,sites,op.first,op.second,{"isHermitian",false}) ;
            auto cm_autompo=AutoMPOCorrelationMatrix<Complex>(psi,sites,op);

            checkDiagonalWithExpect<Complex>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op); 
        } //ops loop
        
    }
    // This type of cross correlation makes a non-symmetric correlation matrix.
    // Let the operator test decide that isHermitian should be false.
    SECTION("correlation_matrix Complex Psi, Sz*Sx, S1/2 No QNs, isHermitian=undef ")
    {
        SiteSet   sites = SpinHalf(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(InitState(sites,"Up"),{"Complex=",true}); 
        
        vector<std::pair<string,string>> ops({{"Sx","Sz"}});

        for (auto op:ops)
        {
            auto cm=correlationMatrixC(psi,sites,op.first,op.second);
            auto cm_autompo=AutoMPOCorrelationMatrix<Complex>(psi,sites,op);

            checkDiagonalWithExpect<Complex>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op); 
        } //ops loop
        
    }
    
    SECTION("correlation_matrix Real Psi, Real Ops, S1/2 With QNs")
    {
        SiteSet   sites = SpinHalf(Nsmall, {"ConserveQNs=",true});
        MPS       psi   = randomMPS(InitState(sites,"Up"),{"Complex=",false}); 
        
        vector<std::pair<string,string>> ops({{"Sz","Sz"},{"S+","S-"},{"S-","S+"}});

        for (auto op:ops)
        {
            auto cm=correlationMatrix(psi,sites,op.first,op.second) ;
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op);

            checkDiagonalWithExpect<Real>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
    }


    SECTION("correlation_matrix Real Psi, Real Ops, Fermions No QNs")
    {
        SiteSet   sites = Fermion(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(sites);
        
        vector<std::pair<string,string>> ops({{"N","N"},{"Cdag","C"},{"Adag","A"},{"C","Cdag"},{"A","Adag"}});
    //   vector<std::pair<string,string>> ops({{"A","Cdag"}}); //Known fail
        for (auto op:ops)
        {
            auto cm     =correlationMatrix(psi,sites,op.first,op.second);
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op);

            checkDiagonalWithExpect<Real>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
      
    }
    SECTION("correlation_matrix Real Psi, Real Ops, Fermions with QNs")
    {
        SiteSet   sites = Fermion(Nsmall, {"ConserveQNs=",true});
        MPS       psi   = randomMPS(InitState(sites,"1"));
        
        vector<std::pair<string,string>> ops({{"N","N"},{"Cdag","C"},{"Adag","A"},{"C","Cdag"},{"A","Adag"}});
    //   vector<std::pair<string,string>> ops({{"A","Cdag"}}); //Known fail
        for (auto op:ops)
        {
            auto cm     =correlationMatrix(psi,sites,op.first,op.second);
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op);

            checkDiagonalWithExpect<Real>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
      
    }

    SECTION("correlation_matrix Real Psi, Real symmetric Ops, Electrons No QNs")
    {
        SiteSet   sites = Electron(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(sites);
        
        vector<std::pair<string,string>> ops({{"Ntot","Ntot"},{"Nup","Nup"},{"Ndn","Ndn"},
                                             {"Cdagup","Cup"},{"Adagup","Aup"},
                                             {"Cdn","Cdagdn"},{"Adn","Adagdn"},
                                             {"Sz","Sz"},{"S+","S-"}
                                             });
        for (auto op:ops)
        {
            auto cm     =correlationMatrix(psi,sites,op.first,op.second);
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op);

            checkDiagonalWithExpect<Real>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
      
    }
    SECTION("correlation_matrix Real Psi, Real non-symmetric Ops, Electrons No QNs")
    {
        SiteSet   sites = Electron(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(sites);
        
        vector<std::pair<string,string>> ops({{"Nup","Ndn"},{"Ntot","Nup"},
                                             {"Sz","Ndn"},{"S+","Nup"}
                                             });
        for (auto op:ops)
        {
            auto cm     =correlationMatrix(psi,sites,op.first,op.second);
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op);

            checkDiagonalWithExpect<Real>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
      
    }
    
    SECTION("correlation_matrix Real Psi, Real Ops, Electrons With QNs")
    {
        SiteSet   sites = Electron(Nsmall, {"ConserveQNs=",true});
        MPS       psi   = randomMPS(InitState(sites,"Up"));
        
        vector<std::pair<string,string>> ops({{"Ntot","Ntot"},{"Nup","Nup"},{"Ndn","Ndn"},
                                             {"Cdagup","Cup"},{"Adagup","Aup"},
                                             {"Cdn","Cdagdn"},{"Adn","Adagdn"},
                                             {"Nup","Ndn"},{"Ntot","Nup"},
                                             {"Sz","Sz"},{"S+","S-"},
                                             });
        for (auto op:ops)
        {
            auto cm     =correlationMatrix(psi,sites,op.first,op.second);
            auto cm_autompo=AutoMPOCorrelationMatrix<Real>(psi,sites,op);

            checkDiagonalWithExpect<Real>(psi,sites,op,cm);
            checkMatrices(cm,cm_autompo,op);
        } //ops loop
      
    }
    SECTION("correlation_matrix expect fail with mixed boson*fermion operator")
    {
        SiteSet   sites = Electron(Nsmall, {"ConserveQNs=",false});
        MPS       psi   = randomMPS(sites);
        
        REQUIRE_THROWS(correlationMatrix(psi,sites,"A","Cdag",{"isHermitian",false}));
      
    }
    

}//TEST_CASE("correlationMatrix")
