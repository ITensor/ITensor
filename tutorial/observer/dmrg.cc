#include "dmrg.h"
#include "sites/spinhalf.h"
#include "sites/spinone.h"
#include "autompo.h"
#include "EntropyObserver.h"

using namespace itensor;

int main()
    {
    //
    // Initialize the sites making up the Hilbert space
    //
    int N = 100;
    //auto sites = SpinHalf(N); //make a chain of N spin 1/2's
    auto sites = SpinOne(N); //make a chain of N spin 1's

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model
    //
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = MPO(ampo);

    // Initalize psi to be a random product
    // MPS on the Hilbert space "sites"
    auto psi = MPS(sites);

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto obs = EntropyObserver(psi);
    auto energy = dmrg(psi,H,sweeps,obs,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);

    //
    // Obtain the energy directly from the MPS
    // and H by computing <psi|H|psi>
    //
    printfln("\nUsing psiHphi = %.10f", psiHphi(psi,H,psi) );

    return 0;
    }
