#include "test.h"
#include "eigensolver.h"
#include "hams/heisenberg.h"
#include "model/spinhalf.h"
#include "localmpo.h"
#include <boost/test/unit_test.hpp>

using namespace std;
using boost::format;

struct EigenSolverDefaults
    {
    EigenSolverDefaults()
        { }
    };

BOOST_FIXTURE_TEST_SUITE(EigenSolverTest,EigenSolverDefaults)

TEST(FourSite)
    {
    //Exact 4 site energy is -1.6160254038 from DMRG

    const int N = 4;
    SpinHalf model(N);
    MPO H = Heisenberg(model);

    InitState initState(N);
    for(int i = 1; i <= N; ++i)
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    MPS psi(model,initState);

    LocalMPO<ITensor> PH(H,2);

    ITensor phip;
    psi.position(2);
    PH.position(2,psi);

    ITensor phi1 = psi.AA(2) * psi.AA(3);

    Eigensolver d(9);
    Real En1 = d.davidson(PH,phi1);
    CHECK_CLOSE(En1,-1.1896926208,1E-4);

    cout << endl << endl;
    /*

    ITensor mpoh = H.AA(2)*H.AA(3);
    ITensor phi2 = psi.AA(2)*psi.AA(3);
    Real En2 = doDavidson(phi2,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=2) = %.20f")%En2 << endl;

    cout << endl << endl;

    cout << endl << "Overlap = " << Dot(phi1,phi2) << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(2,phi1,Fromleft);
    psi.position(3);
    PH.position(3,psi);
    ITensor phi3 = psi.AA(3) * psi.AA(4);
    Real En3 = d.davidson(PH,phi3);
    cout << format("Energy from tensor Davidson (b=3) = %.20f")%En3 << endl;

    cout << endl << endl;

    mpoh = H.AA(3)*H.AA(4);
    ITensor phi3m = psi.AA(3)*psi.AA(4);
    Real En3m = doDavidson(phi3m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=3) = %.20f")%En3m << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(3,phi3,Fromright);
    psi.position(3);
    PH.position(2,psi);
    ITensor phi4 = psi.AA(2) * psi.AA(3);
    Real En4 = d.davidson(PH,phi4);
    cout << format("Energy from tensor Davidson (b=2) = %.20f")%En4 << endl;

    cout << endl << endl;

    mpoh = H.AA(2)*H.AA(3);
    ITensor phi4m = psi.AA(2)*psi.AA(3);
    Real En4m = doDavidson(phi4m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=2) = %.20f")%En4m << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(2,phi4,Fromright);
    psi.position(1);
    PH.position(1,psi);

    //With doDavidson
    mpoh = H.AA(1)*H.AA(2);
    ITensor phi5 = psi.AA(1) * psi.AA(2);
    Real En5 = doDavidson(phi5,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=1) = %.20f")%En5 << endl;

    cout << endl << endl;

    ITensor phi6 = psi.AA(1) * psi.AA(2);
    //Print(phi6);
    //Print(PH.L());
    //Print(PH.R());
    //ITensor AB = phi6 * PH.R();
    //AB *= H.AA(2);
    //AB *= H.AA(1);
    //AB.noprime();
    //ITensor AB; PH.product(phi6,AB);
    //Print(Dot(phi6,AB));
    Real En6 = d.davidson(PH,phi6);
    cout << format("Energy from tensor Davidson (b=1) = %.20f")%En6 << endl;
*/

    }

TEST(IQFourSite)
    {
    //Exact 4 site energy is -1.6160254038 from DMRG

    const int N = 4;
    SpinHalf model(N);
    IQMPO H = Heisenberg(model);

    InitState initState(N);
    for(int i = 1; i <= N; ++i)
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    IQMPS psi(model,initState);

    LocalMPO<IQTensor> PH(H,2);

    IQTensor phip;
    psi.position(2);
    PH.position(2,psi);

    IQTensor phi1 = psi.AA(2) * psi.AA(3);

    Eigensolver d(9);
    Real En1 = d.davidson(PH,phi1);
    cout << format("Energy from tensor Davidson (b=2) = %.20f")%En1 << endl;
    CHECK_CLOSE(En1,-1.1896926208,1E-4);


    }

BOOST_AUTO_TEST_SUITE_END()
