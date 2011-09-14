#define THIS_IS_MAIN
#include "core.h"
#include "hams.h"
#include "dmrg.h"



int main(int argc, char* argv[])
{
    int N = 100;
    int nsweep = 5; 
    int minm = 50;
    int maxm = 100;
    //Real cutoff = 1E-5;
    Real cutoff = 0;

    SpinHalf::Model model(N);

    MPO H = SpinHalf::SquareLattice::Heisenberg(model,1)();

    InitState initState(N);
    for(int i = 1; i <= N; ++i) initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    MPS psi(model,initState);

    cout << boost::format("Initial energy = %.5f\n")%psiHphi(psi,H,psi);

    Sweeps sweeps(ramp_m,nsweep,minm,maxm,cutoff);
    DMRGOpts opts; opts.quiet = true;
    Real En = onesitedmrg(psi,H,sweeps,opts);

    cout << boost::format("\nGround State Energy = %.10f\n")%En;

    return 0;
}
