//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using std::vector;
using std::move;
using namespace itensor;

struct TGate
    {
    int i1 = 0;
    int i2 = 0;
    ITensor G;

    TGate() { }
    TGate(int i1_, int i2_, ITensor G_) 
        : i1(i1_), i2(i2_), G(G_) { }
    };

int
main()
    {
    int N = 20;

    auto sites = SpinHalf(N);

    //Set psi to be Neel state
    auto init = InitState(sites);
    for(auto n : range1(N))
        {
        init.set(n, n%2 == 1 ? "Up" : "Dn");
        }
    auto psi = MPS(init);

    Real ttotal = 2;
    Real tstep = 0.1;

    //Create Trotter gates (imaginary time)
    auto gates = vector<TGate>{};

    for(int b = 1; b < N; ++b)
        {
        auto hh = op(sites,"Sz",b)*op(sites,"Sz",b+1);
        hh += 0.5*op(sites,"Sp",b)*op(sites,"Sm",b+1);
        hh += 0.5*op(sites,"Sm",b)*op(sites,"Sp",b+1);

        auto G = expHermitian(hh,-tstep/2.);

        gates.emplace_back(b,b+1,move(G));
        }
    for(int b = N-1; b >= 1; --b)
        {
        ITensor hh = op(sites,"Sz",b)*op(sites,"Sz",b+1);
        hh += 0.5*op(sites,"Sp",b)*op(sites,"Sm",b+1);
        hh += 0.5*op(sites,"Sm",b)*op(sites,"Sp",b+1);

        auto G = expHermitian(hh,-tstep/2.);

        gates.emplace_back(b,b+1,move(G));
        }

    auto nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(std::fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    for(int step = 1; step <= nt; ++step)
        {
        for(auto& gate : gates)
            {
            auto b = gate.i1;
            auto& G = gate.G;

            psi.position(b);
            auto AA = psi(b)*psi(b+1);

            //
            // TODO: ADD CODE here that applies 
            // the gate G to the MPS bond
            // tensor "AA" by multiplying
            // G and AA using the * operator
            //
            // G is an ITensor
            // with index structure:
            //
            //   s_{b}' s_{b+1}'
            //    |      |
            //    ========
            //    |      |
            //   s_{b}  s_{b+1}
            //
            // After applying G to AA, don't forget
            // to reset the prime level to 0 by using
            // the noPrime method.
            //


            //Normalize AA after applying G
            AA /= norm(AA);

            //SVD AA to restore MPS form
            auto [U,D,V] = svd(AA,inds(psi(b)),{"Cutoff",1E-10});
            psi.set(b,U);
            psi.set(b+1,D*V);
            }

        printfln("Step %d/%d",step,nt);
        }


    //Make Heisenberg H to
    //conveniently measure energy
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = toMPO(ampo);

    printfln("Energy = %.20f",inner(psi,H,psi));

    //
    // Exact ground state energy of N=20
    // Heisenberg model:
    // E0 = -8.6824733306
    //


    return 0;
    }
