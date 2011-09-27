#define THIS_IS_MAIN
#include "core.h"
#include "hams.h"
#include "DMRGWorker.h"

using boost::format;
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char* argv[])
{
    int N = 100;
    int nsweep = 5; 
    int minm = 1;
    int maxm = 100;
    Real cutoff = 1E-5;

    SpinOne::Model model(N);

    MPO H = SpinOne::Heisenberg(model)();

    InitState initState(N);
    for(int i = 1; i <= N; ++i) initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    MPS psi(model,initState);

    cout << format("Initial energy = %.5f\n")%psiHphi(psi,H,psi);

    Sweeps sweeps(Sweeps::ramp_m,nsweep,minm,maxm,cutoff);
    
    DMRGWorker<> dmrgWorker(sweeps);
    
    dmrgWorker.run(psi,H);
    
    Real En = dmrgWorker.energy();
    cout << format("\nGround State Energy = %.10f\n")%En;

    return 0;
}
