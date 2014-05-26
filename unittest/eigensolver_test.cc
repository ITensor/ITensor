#include "test.h"
#include "eigensolver.h"
#include "hams/Heisenberg.h"
#include "sites/spinhalf.h"
#include "localmpo.h"

using namespace itensor;
using namespace std;

TEST_CASE("EigenSolverTest")
{

SECTION("FourSite")
    {
    //Exact 4 site energy is -1.6160254038 from DMRG

    const int N = 4;
    SpinHalf sites(N);
    MPO H = Heisenberg(sites);

    InitState initState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");

    MPS psi(initState);

    LocalMPO<ITensor> PH(H);

    ITensor phip;
    psi.position(2);
    PH.position(2,psi);

    ITensor phi1 = psi.A(2) * psi.A(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    CHECK_CLOSE(En1,-0.95710678118,1E-4);

    cout << endl << endl;
    /*

    ITensor mpoh = H.A(2)*H.A(3);
    ITensor phi2 = psi.A(2)*psi.A(3);
    Real En2 = doDavidson(phi2,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=2) = %.20f")%En2 << endl;

    cout << endl << endl;

    cout << endl << "Overlap = " << Dot(phi1,phi2) << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(2,phi1,Fromleft);
    psi.position(3);
    PH.position(3,psi);
    ITensor phi3 = psi.A(3) * psi.A(4);
    Real En3 = d.davidson(PH,phi3);
    cout << format("Energy from tensor Davidson (b=3) = %.20f")%En3 << endl;

    cout << endl << endl;

    mpoh = H.A(3)*H.A(4);
    ITensor phi3m = psi.A(3)*psi.A(4);
    Real En3m = doDavidson(phi3m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=3) = %.20f")%En3m << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(3,phi3,Fromright);
    psi.position(3);
    PH.position(2,psi);
    ITensor phi4 = psi.A(2) * psi.A(3);
    Real En4 = d.davidson(PH,phi4);
    cout << format("Energy from tensor Davidson (b=2) = %.20f")%En4 << endl;

    cout << endl << endl;

    mpoh = H.A(2)*H.A(3);
    ITensor phi4m = psi.A(2)*psi.A(3);
    Real En4m = doDavidson(phi4m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=2) = %.20f")%En4m << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(2,phi4,Fromright);
    psi.position(1);
    PH.position(1,psi);

    //With doDavidson
    mpoh = H.A(1)*H.A(2);
    ITensor phi5 = psi.A(1) * psi.A(2);
    Real En5 = doDavidson(phi5,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=1) = %.20f")%En5 << endl;

    cout << endl << endl;

    ITensor phi6 = psi.A(1) * psi.A(2);
    //Print(phi6);
    //Print(PH.L());
    //Print(PH.R());
    //ITensor AB = phi6 * PH.R();
    //AB *= H.A(2);
    //AB *= H.A(1);
    //AB.noprime();
    //ITensor AB; PH.product(phi6,AB);
    //Print(Dot(phi6,AB));
    Real En6 = d.davidson(PH,phi6);
    cout << format("Energy from tensor Davidson (b=1) = %.20f")%En6 << endl;
*/

    }

SECTION("IQFourSite")
    {
    //Exact 4 site energy is -1.6160254038 from DMRG

    const int N = 4;
    SpinHalf sites(N);
    IQMPO H = Heisenberg(sites);

    InitState initState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");

    IQMPS psi(initState);

    LocalMPO<IQTensor> PH(H);

    IQTensor phip;
    psi.position(2);
    PH.position(2,psi);

    IQTensor phi1 = psi.A(2) * psi.A(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    //cout << format("Energy from tensor Davidson (b=2) = %.20f")%En1 << endl;
    CHECK_CLOSE(En1,-0.95710678118,1E-4);


    }

}
