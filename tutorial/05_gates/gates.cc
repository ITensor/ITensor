#include "itensor/all.h"
#include "Heisenberg.h"

using std::vector;
using namespace itensor;

int
main()
    {
    int N = 20;

    auto sites = SpinHalf(N);

    auto psi = MPS(sites);

    Real ttotal = 10;
    Real tstep = 0.1;

    auto gates = vector<Gate>();
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
    if(std::fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    for(int step = 1; step <= nt; ++step)
        {
        for(auto& G : gates)
            {
            auto b = G.i1();
            psi.position(b);
            ITensor AA = psi.A(b)*psi.A(b+1);

            //
            // TODO: ADD CODE here that applies 
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

            auto U = psi.A(b);
            ITensor D,V;
            svd(AA,U,D,V,{"Cutoff",1E-10});
            psi.setA(b,U);
            psi.setA(b+1,D*V);
            }
        psi.normalize();
        printfln("Step %d/%d",step,nt);
        }

    MPO H = Heisenberg(sites);
    printfln("Energy = %.20f",overlap(psi,H,psi));


    return 0;
    }
