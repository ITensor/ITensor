#include "tevol.h"
#include "sites/spinhalf.h"
#include "hams/Heisenberg.h"

using std::vector;
using namespace itensor;

int
main(int argc, char* argv[])
    {
    const int N = 20;

    SpinHalf sites(N);

    MPS psi(sites);

    Real ttotal = 10;
    Real tstep = 0.1;

    vector<Gate> gates;
    const Gate::Type type = Gate::tImag;

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

    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    for(int step = 1; step <= nt; ++step)
        {
        Foreach(const Gate& G, gates)
            {
            const int b = G.i();
            psi.position(b);
            const ITensor AA = psi.A(b)*psi.A(b+1);

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



            ITensor D;
            svd(AA,psi.Anc(b),D,psi.Anc(b+1));
            psi.Anc(b+1) *= D;
            }
        psi.normalize();
        printfln("Step %d/%d",step,nt);
        }

    MPO H = Heisenberg(sites);
    printfln("Energy = %.20f",psiHphi(psi,H,psi));


    return 0;
    }
