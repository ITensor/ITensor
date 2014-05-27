#include "sites/spinone.h"
#include "idmrg.h"
#include "hams/Heisenberg.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    //ITensor iDMRG code works with two unit cells
    //in the central region
    const int Nuc = 4;
    const int N = 2*Nuc;

    SpinOne sites(N);

    IQMPO H = Heisenberg(sites,"Infinite=true");

    Sweeps sweeps(20);
    sweeps.maxm() = 20,80,140,200;
    sweeps.cutoff() = 1E-10,Opt("Repeat",10),1E-14;
    sweeps.niter() = 3,2;

    InitState initState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            initState.set(i,"Up");
        else
            initState.set(i,"Dn");
        }
    IQMPS psi(initState);

    idmrgRVal<IQTensor> res = 
    idmrg(psi,H,sweeps,"OutputLevel=2");

    printfln("\nGround state energy / site = %.20f",res.energy/N);

    //Interesting to compare ground state energy to White, Huse, PRB 48, 3844 (1993).

    return 0;
    }
