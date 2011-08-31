#define THIS_IS_MAIN
#include "dmrg.h"
#include "hams.h"

int main(int argc, char* argv[])
{
    int N = 100;
    int nsweep = 5; 
    int minm = 1;
    int maxm = 100;
    Real cutoff = 1E-5;

    SpinOne::Model model(N);

    IQMPO H = SpinOne::Heisenberg(model)();

    InitState initState(N);
    for(int i = 1; i <= N; ++i) initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    IQMPS psi(model,initState);

    cout << format("Initial energy = %.5f\n")%psiHphi(psi,H,psi);

    Sweeps sweeps(ramp_m,nsweep,minm,maxm,cutoff);
    Real En = dmrg(psi,H,sweeps);

    cout << format("\nGround State Energy = %.10f\n")%En;
    cout << "\nTotal QN of Ground State = " << total_QN(psi) << "\n";

    return 0;
}
