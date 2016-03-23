#include "itensor/mps/dmrg.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/mps/autompo.h"

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
    // next-neighbor Heisenberg model.
    //
    // Here we convert the AutoMPO information
    // into an IQMPO, a matrix-product operator
    // which automatically tracks quantum
    // number information.
    //
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = IQMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a classical Neel state
    auto initState = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            initState.set(i,"Up");
        else
            initState.set(i,"Dn");
        }

    //Because we made the initial state an
    //IQMPS with total Sz=0, DMRG will keep
    //it in this same quantum number sector
    auto psi = IQMPS(initState);

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    //Here we use the "noise" feature to aid
    //convergence during the first two sweeps
    sweeps.noise() = 1E-7,1E-8,0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);

    printfln("\nUsing psiHphi = %.10f", psiHphi(psi,H,psi) );

    println("\nTotal QN of Ground State = ",totalQN(psi));

    return 0;
    }
