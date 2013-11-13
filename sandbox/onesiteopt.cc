#define THIS_IS_MAIN
#include "core.h"
#include "model/spinhalf.h"
#include "hams/Heisenberg.h"
using boost::format;
using namespace std;

int main(int argc, char* argv[])
    {
    int N = 100;
    int nsweep = 5; 
    int minm = 50;
    int maxm = 100;
    //Real cutoff = 1E-5;
    Real cutoff = 0;

    SpinHalf model(N);

    MPO H = Heisenberg(model);

    InitState initState(model);
    for(int i = 1; i <= N; ++i) 
        initState.set(i,(i%2==1 ? "Up" : "Dn"));

    MPS psi(initState);

    cout << format("Initial energy = %.5f\n")%psiHphi(psi,H,psi);

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. Here less than 10 maxm
    // values are provided, so all remaining sweeps will use the
    // last maxm (= 200).
    //
    Sweeps sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-8;
    cout << sweeps;

    Real En = onesitedmrg(psi,H,sweeps,"Quiet");

    cout << format("\nGround State Energy = %.10f\n")%En;

    return 0;
    }
