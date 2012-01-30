#define THIS_IS_MAIN
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
    int maxm = 200;
    Real cutoff = 1E-8;

    //SpinHalf model(N);
    SpinOne model(N);

    MPO H = Heisenberg(model);

    InitState initState(N);
    for(int i = 1; i <= N; ++i) 
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    MPS psi(model,initState);

    cout << format("Initial energy = %.5f")%psiHphi(psi,H,psi) << endl;

    Sweeps sweeps(Sweeps::ramp_m,nsweep,minm,maxm,cutoff);

    Real En = dmrg(psi,H,sweeps);

    cout << format("\nGround State Energy = %.10f")%En << endl;

    return 0;
    }
