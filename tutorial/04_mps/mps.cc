//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using std::vector;

int main()
    {
    int N = 50;

    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    //
    // TODO
    //
    // 3. Adjust the external field to see how the
    //    magnetization changes
    //
    Real h = 0.5;

    // Create the MPO for the transverse field
    // Ising model
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += -4.0,"Sz",j,"Sz",j+1;
        }
    for(int j = 1; j <= N; ++j)
        {
        ampo += -2*h,"Sx",j;
        }
    auto H = toMPO(ampo);

    // Create a random starting state
    // For DMRG
    auto psi0 = randomMPS(sites);

    // Run DMRG to get the ground state
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 5,10,20;
    sweeps.cutoff() = 1E-10;
    auto [E,psi] = dmrg(H,psi0,sweeps,{"Quiet",true});
    println("Ground state energy = ",E);

    // Make psidag, the conjugate of psi
    auto psidag = dag(prime(psi));

    // A vector holding the operators used
    // in the expectation value.
    // All set to identity operators to start.
    // Note: this is one-indexed 
    //       (O[n] is the operator on site n)
    auto O = vector<ITensor>(N+1);
    for(auto j : range1(N))
        {
        O[j] = op(sites,"Id",j);
        }

    // Position we will place our operator
    int Npos = N/2;

    //
    // TODO
    //
    // 1. Add an operator to measure the magnetization
    //    in the z direction at Npos.
    //    HINT: use the op(...) function.
    //          It provides spin operators "Sx", "Sy", "Sz",
    //          you may want to scale your operator
    //          to make it a Pauli matrix
    //

    /* Your code here */

    //
    // TODO
    //
    // 2. Complete the following code 
    //    to measure the magnetization.
    //    Print your results with PrintData(...)
    //

    auto o = psidag(1)*O[1]*psi(1);
    for(auto j : range1(2,N))
      {
      /* Your code here */
      }
    
    //
    // TODO
    //
    // 3. Adjust the transverse field h at the top of the
    //    file to find the critical point.
    //    HINT: think about the field limits h -> 0
    //          and h -> infinity
    //

    return 0;
    }
