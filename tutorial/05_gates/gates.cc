#include "tevol.h"
#include "sites/spinhalf.h"
#include "hams/Heisenberg.h"

using std::vector;
using namespace itensor;

int
main(int argc, char* argv[])
    {
    int N = 20;

    SpinHalf sites(N);

    MPS psi(sites);

    Real ttotal = 10;
    Real tstep = 0.1;

    vector<Gate> gates;
    auto type = Gate::tImag;

    for(int b = 1; b < N; ++b)
        {
        ITensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
        hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
        hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
        gates.push_back(Gate(sites,b,b+1,type,tstep/2.,hh));
        }
    for(int b = N-1; b >= 1; --b)
        {
        ITensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
        hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
        hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
        gates.push_back(Gate(sites,b,b+1,type,tstep/2.,hh));
        }

    auto nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    for(int step = 1; step <= nt; ++step)
        {
        for(auto& G : gates)
            {
            auto b = G.i();
            psi.position(b);
            auto AA = psi.A(b)*psi.A(b+1);

            //
            // Write code here that applies 
            // the gate G to the MPS bond
            // tensor "AA"
            //
            // G can be treated as an ITensor
            // with index structure:
            //
            //   s_b'   s_{b+1}'
            //    |      |
            //    ========
            //    |      |
            //   s_b    s_{b+1}
            //
            // After applying G to AA, don't forget
            // to reset the prime level to 0 by using
            // the noprime or mapprime methods.
            //



            ITensor D;
            svd(AA,psi.Anc(b),D,psi.Anc(b+1));
            psi.Anc(b+1) *= D;
            }
        psi.normalize();
        printfln("Step %d/%d",step,nt);
        }

    auto H = MPO(Heisenberg(sites));
    printfln("Energy = %.20f",psiHphi(psi,H,psi));


    return 0;
    }
