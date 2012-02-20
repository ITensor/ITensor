#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/heisenberg.h"
using boost::format;
using namespace std;

int main(int argc, char* argv[])
    {
    int N = 100;
    int nsweep = 5;
    int minm = 1;
    int maxm = 100;
    Real cutoff = 1E-8;

    //SpinHalf model(N);
    SpinOne model(N);

    IQMPO H = Heisenberg(model);

    InitState initState(N);
    for(int i = 1; i <= N; ++i) 
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    IQMPS psi(model,initState);

    cout << format("Initial energy = %.5f\n")%psiHphi(psi,H,psi);

    Sweeps sweeps(Sweeps::ramp_m,nsweep,minm,maxm,cutoff);

    Real En = dmrg(psi,H,sweeps);

    cout << format("\nGround State Energy = %.10f\n")%En;
    cout << "\nTotal QN of Ground State = " << totalQN(psi) << "\n";

    return 0;
    }
